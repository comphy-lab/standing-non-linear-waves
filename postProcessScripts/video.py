# Author: Vatsal Sanjay
# vatsalsanjay@gmail.com
# Physics of Fluids
# Last updated: Dec 24, 2024

import numpy as np
import os
import subprocess as sp
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.ticker import StrMethodFormatter
import multiprocessing as mp
from functools import partial
import argparse  # Add at top with other imports

import matplotlib.colors as mcolors
custom_colors = ["white", "#DA8A67", "#A0522D", "#400000"]
custom_cmap = mcolors.LinearSegmentedColormap.from_list("custom_hot", custom_colors)

matplotlib.rcParams['font.family'] = 'serif'
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

def gettingFacets(filename):
    exe = ["./getFacets", filename]
    p = sp.Popen(exe, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    temp1 = stderr.decode("utf-8")
    temp2 = temp1.split("\n")
    segs = []
    skip = False
    if (len(temp2) > 1e2):
        for n1 in range(len(temp2)):
            temp3 = temp2[n1].split(" ")
            if temp3 == ['']:
                skip = False
                pass
            else:
                if not skip:
                    temp4 = temp2[n1+1].split(" ")
                    r1, z1 = np.array([float(temp3[1]), float(temp3[0])])
                    r2, z2 = np.array([float(temp4[1]), float(temp4[0])])
                    segs.append(((r1, z1),(r2, z2)))
                    skip = True
    return segs

def gettingfield(filename, zmin, rmin, zmax, rmax, nr):
    exe = ["./getData", filename, str(zmin), str(rmin), str(zmax), str(rmax), str(nr)]
    p = sp.Popen(exe, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = p.communicate()
    temp1 = stderr.decode("utf-8")
    temp2 = temp1.split("\n")
    # print(temp2) #debugging
    Rtemp, Ztemp, D2temp, veltemp, omegaTemp, fTemp, uxTemp, uyTemp  = [],[],[],[],[],[],[],[]

    for n1 in range(len(temp2)):
        temp3 = temp2[n1].split(" ")
        if temp3 == ['']:
            pass
        else:
            Ztemp.append(float(temp3[0]))
            Rtemp.append(float(temp3[1]))
            D2temp.append(float(temp3[2]))
            veltemp.append(float(temp3[3]))
            omegaTemp.append(float(temp3[4]))
            fTemp.append(float(temp3[5]))
            uxTemp.append(float(temp3[6]))
            uyTemp.append(float(temp3[7]))


    R = np.asarray(Rtemp)
    Z = np.asarray(Ztemp)
    D2 = np.asarray(D2temp)
    vel = np.asarray(veltemp)
    omega = np.asarray(omegaTemp)
    f = np.asarray(fTemp)
    ux = np.asarray(uxTemp)
    uy = np.asarray(uyTemp)

    nz = int(len(Z)/nr)

    # print("nr is %d %d" % (nr, len(R))) # debugging
    print("nz is %d" % nz)

    R.resize((nz, nr))
    Z.resize((nz, nr))
    D2.resize((nz, nr))
    vel.resize((nz, nr))
    omega.resize((nz, nr))
    f.resize((nz, nr))
    ux.resize((nz, nr))
    uy.resize((nz, nr))

    return R, Z, D2, vel, omega, f, ux, uy, nz

def save_datasets(filename, X, Y, f_rot, vel_rot, omega_rot, ux_rot, uy_rot):
    """
    Save all relevant datasets into a single .npz file.
    You can read them later with:
      data = np.load('filename.npz')
      X = data['X']
      Y = data['Y']
      ...
    """
    np.savez(filename,
             X=X, Y=Y,
             f_rot=f_rot,
             vel_rot=vel_rot,
             omega_rot=omega_rot,
             ux_rot=ux_rot,
             uy_rot=uy_rot)
# ----------------------------------------------------------------------------------------------------------------------

def process_timestep(ti, caseToProcess, folder, nGFS, GridsPerR, rmin, rmax, zmin, zmax, lw):
    t = 0.01 * ti
    place = f"{caseToProcess}/intermediate/snapshot-{t:.4f}"
    name = f"{folder}/{int(t*1000):08d}.png"

    if not os.path.exists(place):
        print(f"{place} File not found!")
        return

    if os.path.exists(name):
        print(f"{name} Image present!")
        return

    # --- Gather interface segments and field data as before ---
    segs = gettingFacets(place)
    nr = int(GridsPerR * rmax)
    R, Z, D2, vel, omega, f, ux, uy, nz = gettingfield(place, zmin, rmin, zmax, rmax, nr)

    # ----------------------------------------------------------
    # #FIXME: #TODO: this is not elegant. Ensure that the variables are all taken properly from the very beginning... 
    # ----------------------------------------------------------
    # 1) Prepare arrays for a 90° rotation: we want x=Z, y=R.
    #    Currently, vel.shape == (nz, nr) with vel[iZ, iR].
    #    After transposing, shape becomes (nr, nz), so index 0
    #    will correspond to R and index 1 to Z. We also do it
    #    for omega:
    # ----------------------------------------------------------
    X = Z.T
    Y = R.T
    vel_rot   = vel.T
    omega_rot = omega.T
    f_rot     = f.T
    ux_rot    = ux.T
    uy_rot    = uy.T

    # Save data to disk
    out_datafile = f"{folder}/{int(t*1000):08d}_data.npz"
    save_datasets(out_datafile, X, Y, f_rot, vel_rot, omega_rot, ux_rot, uy_rot)

    # ----------------------------------------------------------
    # 2) Swap the meaning of x<->z and y<->r in the extent:
    #    old extent was [rmin, rmax, zmin, zmax].
    #    now we want x from zmin to zmax, y from rmin to rmax.
    # ----------------------------------------------------------
    x_extent = [zmin, zmax]
    y_extent = [rmin, rmax]

    # For imshow, extent is [x_min, x_max, y_min, y_max].
    extent_vel   = [x_extent[0], x_extent[1], y_extent[0], y_extent[1]]
    extent_omega = [x_extent[0], x_extent[1], y_extent[0], y_extent[1]]

    # ----------------------------------------------------------
    # 3) Rotate the interface segments: each segment is ((r1,z1),(r2,z2)).
    #    We want ((z1,r1),(z2,r2)) to swap x<->z, y<->r.
    # ----------------------------------------------------------
    segs_rot = []
    for seg in segs:
        (r1, z1), (r2, z2) = seg
        segs_rot.append(((z1, r1), (z2, r2)))

    # ----------------------------------------------------------
    # 4) Make subplots.  Often, one does subplots(1,2) if one
    #    wants side-by-side panels. If you prefer them stacked,
    #    keep subplots(2,1).  Either is fine, as long as you
    #    realize you have “rotated” the data. For side-by-side:
    # ----------------------------------------------------------
    AxesLabel, TickLabel = 50, 20
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(19.20, 10.80))

    # Draw the rotated interface in green, domain boundaries, etc.
    for ax in [ax1, ax2]:
        # Gray dashed lines now get swapped similarly if you want
        # them at z=0 or r=0.  If you simply want a bounding box,
        # note that x->z, y->r:
        ax.plot([zmin, zmax], [0, 0], '--', color='grey', linewidth=lw)  # "horizontal" axis is z
        ax.plot([0, 0], [rmin, rmax], '-.', color='grey', linewidth=lw)  # "vertical" axis is r

        # Domain box:
        ax.plot([zmin, zmax], [rmin, rmin], '-', color='black', linewidth=lw)
        ax.plot([zmin, zmax], [rmax, rmax], '-', color='black', linewidth=lw)
        ax.plot([zmin, zmin], [rmin, rmax], '-', color='black', linewidth=lw)
        ax.plot([zmax, zmax], [rmin, rmax], '-', color='black', linewidth=lw)

        line_segments = LineCollection(segs_rot, linewidths=4, colors='green')
        ax.add_collection(line_segments)

    # ----------------------------------------------------------
    # 5) Now show imshow with the rotated arrays and extents:
    # ----------------------------------------------------------
    cntrl1 = ax1.imshow(
        vel_rot, 
        cmap="Blues", 
        interpolation='bilinear', 
        origin='lower', 
        extent=extent_vel,
        vmin=0.0,
        vmax=0.5
    )
    cntrl2 = ax2.imshow(
        omega_rot, 
        cmap="coolwarm", 
        interpolation='bilinear', 
        origin='lower', 
        extent=extent_omega,
        vmin=-1e0,
        vmax=1e0
    )

    # Equal aspect ensures squares in the new orientation
    for ax in [ax1, ax2]:
        ax.set_aspect('equal')
        ax.set_xlim(zmin, zmax)  # x range
        ax.set_ylim(rmin, rmax)  # y range

    # Titles and labels that match the new orientation
    ax1.set_title(fr'$t\sqrt{{g/\lambda}} = {t:.4f}$ (Velocity)', fontsize=TickLabel)
    # ax1.set_xlabel('$X$', fontsize=AxesLabel)
    # ax1.set_ylabel('$r$', fontsize=AxesLabel)

    ax2.set_title(r'Vorticity', fontsize=TickLabel)
    # ax2.set_xlabel('$z$', fontsize=AxesLabel)
    # ax2.set_ylabel('$r$', fontsize=AxesLabel)

    # Colorbars: place them below each subplot, for instance
    fig.subplots_adjust(bottom=0.2, wspace=0.3)  # more spacing for colorbars
    cbar_ax1 = fig.add_axes([0.125, 0.1, 0.35, 0.03])   # x,y,width,height in figure coords
    c1 = plt.colorbar(cntrl1, cax=cbar_ax1, orientation='horizontal')
    c1.ax.tick_params(labelsize=TickLabel)
    c1.set_label(r'$\|u/\sqrt{g\lambda}\|$', fontsize=AxesLabel)

    cbar_ax2 = fig.add_axes([0.57, 0.1, 0.35, 0.03])
    c2 = plt.colorbar(cntrl2, cax=cbar_ax2, orientation='horizontal')
    c2.ax.tick_params(labelsize=TickLabel)
    c2.set_label(r'$\omega\sqrt{\lambda/g}$', fontsize=AxesLabel)

    for ax in [ax1, ax2]:
        ax.axis('off')

    plt.savefig(name, bbox_inches="tight")
    plt.close()

def main():
    # Get number of CPUs from command line argument, or use all available
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--CPUs', type=int, default=mp.cpu_count(), help='Number of CPUs to use')
    parser.add_argument('--nGFS', type=int, default=1550, help='Number of restart files to process')

    parser.add_argument('--ZMAX', type=float, default=1.0, help='Maximum Z value')
    parser.add_argument('--RMAX', type=float, default=1.0, help='Maximum R value')
    parser.add_argument('--ZMIN', type=float, default=-1.0, help='Minimum Z value')
    parser.add_argument('--RMIN', type=float, default=-1.0, help='Minimum R value')
    parser.add_argument('--caseToProcess', type=str, default='StokesStandingWaves', help='Case to process')
    
    args = parser.parse_args()

    CPUStoUse = args.CPUs
    nGFS = args.nGFS
    ZMAX = args.ZMAX
    RMAX = args.RMAX
    ZMIN = args.ZMIN
    RMIN = args.RMIN
    caseToProcess = args.caseToProcess
    num_processes = CPUStoUse
    
    rmin, rmax, zmin, zmax = [RMIN, RMAX, ZMIN, ZMAX]
    GridsPerR = 128

    lw = 2
    folder = 'Video'

    if not os.path.isdir(folder):
        os.makedirs(folder)

    # Create a pool of worker processes
    with mp.Pool(processes=num_processes) as pool:
        # Create partial function with fixed arguments
        process_func = partial(process_timestep, caseToProcess=caseToProcess, 
                             folder=folder, nGFS=nGFS,
                             GridsPerR=GridsPerR, rmin=rmin, rmax=rmax, 
                             zmin=zmin, zmax=zmax, lw=lw)
        # Map the process_func to all timesteps
        pool.map(process_func, range(nGFS))

if __name__ == "__main__":
    main()
