# randomizer.py
# Jasper Lankhorst's thesis
# Reads, centers and randomizes SiQDs

import random
import math
import numpy as np


def reader(input_file):
    """
    Reads inputfile, makes dict of all coordinates. input_file has to be a VESTA
    .xyz file. First two lines are skipped. Return list of dictionaries with
    coordinates.
    """
    lines = open(input_file, "r")

    # skip first two lines
    next(lines)
    next(lines)

    # read lines, add to list of coordinates
    atomlist = []
    for i in lines:
        element, u, v, w = i.split()
        u = float(u)
        v = float(v)
        w = float(w)
        coor_dict = {"name": element, "x": u, "y": v, "z": w}
        atomlist.append(coor_dict)

    return atomlist


def shifter(list_of_atoms):
    """
    Takes list with dicts of coordinates, centers it around (0,0,0).
    Returns the approximate radius.
    """

    # find extremes
    x_max, x_min, y_max, y_min, z_max, z_min = 0, 0, 0, 0, 0, 0
    for atom in list_of_atoms:
        x, y, z  = atom["x"], atom["y"], atom["z"]

        if x > x_max:
            x_max = x
        if x < x_min:
            x_min = x
        if y > y_max:
            y_max = y
        if y < y_min:
            y_min = y
        if z > z_max:
            z_max = z
        if z < z_min:
            z_min = z

    # find offset
    x_offset = x_max + x_min
    y_offset = y_max + y_min
    z_offset = z_max + z_min

    # shift every atom by the negative of the offset, centering it
    for atom in list_of_atoms:
        atom["x"] -= x_offset
        atom["y"] -= y_offset
        atom["z"] -= z_offset

    return max([x_max, abs(x_min), y_max, abs(y_min), z_max, abs(z_min)])


def renamer(atomlist, size, d):
    """
    Rename the atoms in the shell to "Si2". d is the fraction of the radius of
    the QD that is shell.
    """
    amount_inside = 0
    core_r = (size / 2) * (1 - d)
    for atom in atomlist:
        x,y,z = atom["x"], atom["y"], atom["z"]
        if math.sqrt(x**2 + y**2 + z**2) > core_r and atom["name"] == "Si":
            atom["name"] = "Si2"

    return atomlist


def edge_finder(atomlist, size):
    """
    Finds all atoms on the edge of the QD. Renames to "Si3".
    """
    # Make variables for the corners
    # Also: mark every point in the edge you find as Si3
    x_corner_y, x_corner_z, y_corner_x, z_corner_x, z_corner_y, \
    y_corner_z = 0, 0, 0, 0, 0, 0
    for atom in atomlist:
        if abs(atom["x"]) == size/2:
            atom["name"] = "Si3"
            if abs(atom["y"]) >= y_corner_x:
                y_corner_x = abs(atom["y"])
            if abs(atom["z"]) >= z_corner_x:
                z_corner_x = abs(atom["z"])
        if abs(atom["y"]) == size/2:
            atom["name"] = "Si3"
            if abs(atom["x"]) >= x_corner_y:
                x_corner_y = abs(atom["x"])
            if abs(atom["z"]) >= z_corner_y:
                z_corner_y = abs(atom["z"])
        if abs(atom["z"]) == size/2:
            atom["name"] = "Si3"
            if abs(atom["y"]) >= y_corner_z:
                y_corner_z = abs(atom["y"])
            if abs(atom["x"]) >= x_corner_z:
                x_corner_z = abs(atom["x"])

    # Now we know all coordinates of the corners, since the QD is symmetric in
    # all directions.

    # Coordinates of middle point of the normal vector of the 'triangle'
    nvc = (size/2 + y_corner_x + z_corner_x) / 3
    normal_vector = [nvc, nvc, nvc]

    # Mark remaining atoms
    for atom in atomlist:
        if atom["name"] == "Si" or atom["name"] == "Si2":
            x, y, z = atom["x"], atom["y"], atom["z"]
            plane_factor = (size - x) + (y_corner_x - y) + (z_corner_x - z)
            if plane_factor < 0.05:
                atom["name"] = "Si3"


def vector(max_length):
    """
    Generates a random direction max_length vector in 3D. Output: x, y, z.
    """
    theta = random.random() * math.pi
    phi = random.random() * 2 * math.pi
    r = max_length
    x = r * math.sin(theta) * math.cos(phi)
    y = r * math.sin(theta) * math.sin(phi)
    z = r * math.cos(theta)
    return x, y, z


def randomizer(atomlist, shift, size, nominal_size, fraction):
    """
    Takes list of dicts of coordinates, shifts them all a length of shift (in
    angstrom) times a random . Every atom's distance from core is multiplied
    by 'fraction'.
    """
    for atom in atomlist:
        if atom["name"] == "Si2" or atom["name"] == "Si3":
            x,y,z = atom["x"], atom["y"], atom["z"]
            dx, dy, dz = vector(np.random.normal(loc=shift,scale=shift/12.0))
            r_old = math.sqrt(x**2 + y**2 + z**2)
            x += dx
            y += dy
            z += dz
            r_new = math.sqrt(x**2 + y**2 + z**2)
            if r_new > r_old:
                x -= 2 * dx
                y -= 2 * dy
                z -= 2 * dz

            atom["x"], atom["y"], atom["z"] = pusher_backer(x, y, z, fraction)

    return atomlist


def saver(atomlist, percentage, size, d, filename):
    """
    Saves the atomlist in a new .xyz, n is the iteration number.
    """
    # shell = 1 - d
    amount = len(atomlist)
    tempfile = '/home/jasper/Desktop/studie/thesis/Code/Si{}_shell{}R_{}%.xyz'\
    .format(filename, d, percentage)
    with open(tempfile, "w+") as file:
        file.write("{}\n \n".format(str(amount)))
        for atom in atomlist:
            name, x, y, z = atom["name"], atom["x"], atom["y"], atom["z"]
            atomstring = "{}  {}  {}  {} \n".format(name, x, y, z)
            file.write(atomstring)


def pusher_backer(x, y, z, fraction):
    """
    Distance from core times 'fraction'.
    """
    x *= fraction
    y *= fraction
    z *= fraction
    return x, y, z


def main(filepath, d, bond_length, percentage, iterations, filename,
            nominal_size, fraction):
    """
    Main function. Saves in new .xyz.
    """
    lijst = reader(filepath)
    size = 2 * shifter(lijst)
    shift = percentage / 100 * bond_length

    renamer(lijst, size, d)

    edge_finder(lijst, size)

    if d > 0:
        lijst = randomizer(lijst, shift, size, nominal_size, fraction)

    saver(lijst, percentage, size, d, filename)


if __name__ == '__main__':
    filerange = [29, 227, 453, 1947]
    d_range = [0, 0.1, 0.4, 1] # fraction of radius that is shell (not core)
    bond_length = 2.361
    percentage = 35
    iterations = 1
    sizerange = [1, 2, 3, 4]
    fraction = 0.9
    for d in d_range:
        for nominal_size in sizerange:
            filename = filerange[int(nominal_size) - 1]
            main("/home/jasper/Desktop/studie/thesis/VESTA-gtk3/Si{}.xyz"\
            .format(filename), d, bond_length, percentage, iterations,\
            filename, nominal_size, fraction)
