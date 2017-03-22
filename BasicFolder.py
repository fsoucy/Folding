import pyrosetta
import pyrosetta.teaching
import rosetta.core.fragment as fragment
import rosetta.protocols.simple_moves as simple_moves
import numpy as np
import time
import random
import math
pyrosetta.init()


def fragment_mover_init(frag_file, frag_length=3):
    fragset = fragment.ConstantLengthFragSet(frag_length)
    fragset.read_fragment_file(frag_file)

    movemap = pyrosetta.MoveMap()
    movemap.set_bb(True)
    fragment_mover = simple_moves.ClassicFragmentMover(fragset, movemap)

    return fragment_mover


def single_move(pose):
    res = random.randint(a=1, b=len(pose.sequence()))
    pick_angle = np.random.randint(low=0, high=1, size=1)
    if pick_angle:
        new_angle = random.gauss(mu=pose.phi(res), sigma=25)
        pose.set_phi(res, new_angle)
    else:
        new_angle = random.gauss(mu=pose.psi(res), sigma=25)
        pose.set_psi(res, new_angle)

def angles_from_file(three_mer_file):
    positions = {}
    f = open(three_mer_file, "r")
    lines = f.readlines()
    f.close()
    current = ""
    count = 0
    angles = []
    for l in lines:
        if "position" in l:
            spl = l.split(" ")
            spl = list(filter(lambda l: len(l) > 0, spl))
            pos = int(spl[1])
            current = pos
            positions[current] = []
        else:
            if (not l == '\n'):
                phi, psi, omega = get_angles(l)
                angles.append((phi, psi, omega, l))
                count += 1
                if count == 3:
                    positions[current].append(angles)
                    angles = []
                    count = 0
    return positions

def get_angles(line):
    x = line.split(" ")
    x = list(filter(lambda l: len(l) > 0, x))
    print(x)
    phi, psi, omega = x[5], x[6], x[7][:-1]
    print(phi)
    print(psi)
    return (float(phi), float(psi), float(omega))

def single_frag_move(pose, positions):
    position = random.randint(1, len(pose.sequence()))
    fragments = positions[position]
    fragment = fragments[random.randint(0, len(fragments) - 1)]
    offset = 0
    frag_description = "\n"
    while offset < 3:
        phi, psi, omega, description = fragment[offset]
        frag_description += description + "\n"
        pose.set_phi(position + offset, phi)
        pose.set_psi(position + offset, psi)
        pose.set_omega(position + offset, omega)
        offset += 1
    return frag_description

seq = 'MIKVTVTNSFFEVTGHAPDKTLCASVSLLTQHVANFLKAEKKAKIKKESGYLKVKFEELENCEVKVLAAMVRSLKELEQKFPSQIRVEVIDNGS'
three_mer_file = '../Projects/homework/homework_folding/structure_prediction/1S12/1S12.200.3mers'
nine_mer_file = '../Projects/homework/homework_folding/structure_prediction/1S12/1S12.200.9mers'
fragments = angles_from_file(three_mer_file)

pose = pyrosetta.pose_from_sequence(seq)
pmm = pyrosetta.PyMOLMover()
pmm.apply(pose)

accepted = 0
rejected = 0

scorefxn = pyrosetta.teaching.get_fa_scorefxn()
score = scorefxn(pose)

best_score = score
best_pose = pyrosetta.pose_from_sequence(seq)
best_pose.assign(pose)

new_pose = pyrosetta.pose_from_sequence(seq)

frag_mover = fragment_mover_init(frag_file=three_mer_file, frag_length=3)

accepted_fragment_file = open("accepted_fragments.txt", "a")
rejected_fragment_file = open("rejected_fragments.txt", "a")

for i in range(10):
    new_pose.assign(pose)
    #frag_mover.apply(new_pose)
    fragment = single_frag_move(new_pose, fragments)
    new_score = scorefxn(new_pose)
    r = random.random()

    if new_score <= score or r < math.exp(score - new_score):
        pose.assign(new_pose)
        new_pose.dump_pdb("accepted_" + str(accepted) + ".pdb")
        accepted_fragment_file.write(fragment)
        score = new_score
        pmm.apply(pose)
        accepted += 1
    # time.sleep(.5)
    else:
        rejected_fragment_file.write(fragment)
        new_pose.dump_pdb("rejected_" + str(rejected) + ".pdb")
        rejected += 1
    if new_score <= best_score:
        best_pose.assign(new_pose)
        best_score = new_score

best_pose.dump_pdb("best_pose.pdb")
