# foldon
Scripts specific to the foldon project that you seriously freaking need to keep track of. 

addloop_renum.py
    adds a loop defined by 'loop.pdb' to an existing pdb file. All atom indices are hard-coded for 1rfo.pdb

alex rmsd finder
    copied from Alex DeGrage (ajd98)

conc.py
    calculates the concentration of your (currently dimer) in a 125A sphere.

count_events_200
    calculates kinetic properties (k_on, off, etc) based on number of events/total time in a state

count_native_contacts
    counts the number of unique native contacts in a given system.go.contacts format.

dihedrals.py
    functional script which modifies internal.parameters based on atom indices

dihedrals_ali
    presumably functional script which modifies many things in 2 chains of a single pdb. written by ali sinan saglam (asinansaglam)

indi comp plots
    plotting script for comparison 

native_contacts
    more or less a copy of ugh_fraction_native_contacts, written with Alex DeGrage (ajd98)

new go contacts
    i literally can't remember what this does

nonrandom
    modified randomorientations (from ajd98) which separates chains of 1rfo and separates them in native orientation

one rmsd
    modified rmsd_2 for a single calculation

plot fr nc
    srsly old plotting script for native contacts array output from ugh_fraction_native_contacts.py

plot moving average
    example plotting script for a moving instantaneous average

randomorientations
    script by ajd98 for orienting monomers according to a random matrix. sources within code because alex is the best

randomorientations restart
    modified version of the above which works with the restart.file format rather than .pdb

renum go contacts ASS
    modified from ali sinan saglam (asinansaglam) to renumber go contacts within system.go.contacts

renum oriented
    takes 1rfo_oriented (output of separatechains) and removes junk. **poorly named**

rmsd
    from ajd98

rmsd_2
    modified from ajd98, calculates rmsd relative to permutations of the three monomers

rmsd scf
    modified of rmsd_2, calculates rmsd relative to reference scf structure

separatechains
    takes 1rfo.pdb and writes 3 pdbs, one for each monomer

tails.py
    wow this was a long time ago. theoretically takes rmsd of tail regions of 1rfo monomers, tracks along course of a trajectory.

tryp_calc
    calculates minimum distance between tryptophan residues in 3 monomers of 1rfo

tryptophan tracker & tryptophan tracker scf
    uses same method as tails.py to track tryptophans through course of trajectory (untested)

ugh fraction native contacts
    calculates and records atom pairs which are within 5.5A +- 1.2*equilibrium distance in a given structure

ugh scf nc
    does the same thing as above, with the atom indices changed for the single chain foldon. 
