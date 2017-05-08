import sys
from rdkit import Chem
from rdkit.Chem import AllChem

max_charge = 2

filename = sys.argv[1]
file = open(filename, "r")

smartsref = ( ('[NX3;H2;!$(NC=O)&!$([NX3+])&!$(NC=N)]','[NH3+]'),
              ('[NX3;H1;!$(NC=O)&!$([NX3+])&!$(NC=N)]','[NH2+]'),
              ('[NX3;H0;!$(NC=O)&!$([NX3+])&!$(NC=N)]','[NH+]'),
#              ('[NX3;H0;!$(NC=O)&!$([NX3+])]','[NH+]'),
              ('[NX2;H1]','[NH2+]'), 
              ('[NX2;H0]','[NH+]'), 
              ('[nX2;H0]','[NH+]'), 
              ('[OX2H1;!$(O[CX4])]','[O-]'),
              ('[SX2H1;!$(S[CX4])]','[S-]') )

def tautomerize(m):
    tautomers = []
    rxn_smarts = ['[C:1]=[NX2;H0:2]-[N;H1:3]-[C:4]>>[C:1]-[N;H1:2]-[N;H0:3]=[C:4]',
                  '[nX2;H0:1]:[n;H1:2]>>[n;H1:1]:[n;H0:2]',
                  '[NX2;H0:1]=[C,N:2]-[N;H1:3]>>[N;H1:1]-[*:2]=[N;H0:3]',
                  '[nX2;H0:1]:[c,n:2]:[n;H1:3]>>[n;H1:1]:[*:2]:[n;H0:3]',
                  '[NX2;H0:1]=[CR0,NR0:2]-[CR0,NR0:3]=[CR0,NR0:4]-[N;H1:5]>>[N;H1:1]-[*:2]=[*:3]-[*:4]=[N;H0:5]',
                  '[NX2;H0:1]=[C,N:2]-[c,n:3]:[c,n:4]-[N;H1:5]>>[N;H1:1]-[*:2]=[*:3]:[*:4]=[N;H0:5]',
                  '[NX2;H0:1]=[c,n:2]:[c,n:3]=[C,N:4]-[N;H1:5]>>[N;H1:1]-[*:2]:[*:3]-[*:4]=[N;H0:5]',
#                 '[nr5;H1:1]:[c,n:2]:[c,n:3]:[c,n:4]:[nr6;H0:5]>>[nr5;H0:1]:[*:2]:[*:3]-[*:4]:[nr6;H1:5]',
                  '[NX2;H0:1]=[C,N:2]-[N;H2:3]>>[N;H1:1]-[*:2]=[N;H1:3]',
                  '[NX2;H1:1]=[C,N:2]-[N;H1:3]>>[N;H2:1]-[*:2]=[N;H0:3]',
#should acetyl be included? What about esters? thione/thiol?
                  '[C;H2:1]=[C:2]-[O;H1:3]>>[C;H3:1]-[*:2]=[O;H0:3]',
                  '[C;H1:1]=[C:2]-[O;H1:3]>>[C;H2:1]-[*:2]=[O;H0:3]',
                  '[C;H0:1]=[C:2]-[O;H1:3]>>[C;H1:1]-[*:2]=[O;H0:3]']
#                 '[C;H3:1]-[C:2](=[O;H0:3])-[C,c:4]>>[C;H1:1]=[*:2](-[O;H1:3])-[C:4]',
#                 '[C;H2:1]-[C:2](=[O;H0:3])-[C,c:4]>>[C;H1:1]=[*:2](-[O;H1:3])-[C:4]',
#                 '[C;H1:1]-[C:2](=[O;H0:3])-[C,c:4]>>[C;H0:1]=[*:2](-[O;H1:3])-[C:4]']
#                  '[N;H0:1]=[C,N:2]-[C,N:3]=[C,N:4]-[C,N:5]=[C,N:6]-[N;H1:7]>>[N;H1:1]-[*:2]=[*:3]-[*:4]=[*:5]-[*:6]=[N;H0:7]',
#                  '[N;H0:1]=[C,N:2]-[C,N:3]=[C,N:4]-[C,N:5]=[C,N:6]-[C,N:7]=[C,N:8]-[N;H1:9]>>[N;H1:1]-[*:2]=[*:3]-[*:4]=[*:5]-[*:6]=[*:7]-[*:8]=[N;H0:9]']
    for smarts in rxn_smarts:
        rxn = AllChem.ReactionFromSmarts(smarts)
        ps = rxn.RunReactants((m,))
        for x in ps:
            smiles = Chem.MolToSmiles(x[0],isomericSmiles=True)
#           print Chem.MolToSmiles(m,isomericSmiles=True),smarts,smiles
            if smiles not in tautomers: tautomers.append(smiles)
     
    for i in range(len(tautomers)):
         for tautomer in tautomers:
             m = Chem.MolFromSmiles(tautomer)
             for smarts in rxn_smarts:
#                print tautomer,Chem.MolToSmiles(m,isomericSmiles=True),smarts
                 rxn = AllChem.ReactionFromSmarts(smarts)
                 ps = rxn.RunReactants((m,))
                 for x in ps:
                     smiles = Chem.MolToSmiles(x[0],isomericSmiles=True)
#                    print tautomer,Chem.MolToSmiles(m,isomericSmiles=True),smarts,smiles
                     if smiles not in tautomers: tautomers.append(smiles)

    return tautomers

def enol_to_keto(smiles):
#
# assumes one enol. Fix for multiple later
#
    m = Chem.MolFromSmiles(smiles)
    rxn_smarts = ['[C;H2:1]=[C:2]-[O;H1:3]>>[C;H3:1]-[*:2]=[O;H0:3]',
                  '[C;H1:1]=[C:2]-[O;H1:3]>>[C;H2:1]-[*:2]=[O;H0:3]',
                  '[C;H0:1]=[C:2]-[O;H1:3]>>[C;H1:1]-[*:2]=[O;H0:3]']

    for smarts in rxn_smarts:
        rxn = AllChem.ReactionFromSmarts(smarts)
        ps = rxn.RunReactants((m,))
        if len(ps) > 0: smiles = Chem.MolToSmiles(ps[0][0],isomericSmiles=True)
        
    return smiles


def in_diff_rings(m,charged_atoms,rings):
    rings=[set(x) for x in rings]
    all_rings = []
    for i in range(len(rings)):
        all_rings.append(rings[i])
        for j in range(i+1,len(rings)):
            overlap=rings[i]&rings[j]
            if len(overlap) > 1:
                all_rings.append(rings[i]|rings[j])

#   print all_rings, charged_atoms
    pair = []
    for ring in all_rings:
        for i in charged_atoms:
            atom_i = int(i[0])
            pair.append(atom_i)
            for j in charged_atoms:
                if j <= i: continue
                atom_j = int(j[0])
                pair.append(atom_j)
#               print pair, list(ring),all(x in list(ring) for x in pair)
                if all(x in list(ring) for x in pair) and m.GetAtomWithIdx(atom_i).GetIsAromatic() and m.GetAtomWithIdx(atom_j).GetIsAromatic():
                   return False
                    
    return True

for line in file:
    ionsmiles = []
    words = line.split()
    name = words[0]
    neutral_smiles = words[1]
#
#canonicalize
#
    neutral_smiles = Chem.MolToSmiles(Chem.MolFromSmiles(neutral_smiles),isomericSmiles=True)
#
# enol -> keto
#
    neutral_smiles = enol_to_keto(neutral_smiles)
#
    ionsmiles.append(neutral_smiles)
    neutral_m = Chem.MolFromSmiles(neutral_smiles)
    tautomers = tautomerize(neutral_m)
#   tautomers = []
#   print tautomers
   
    tautomers.insert(0, neutral_smiles)
 
    for tautomer in tautomers:
        if tautomer not in ionsmiles: ionsmiles.append(tautomer)
        m = Chem.MolFromSmiles(tautomer)

        i = 0
        for (smarts1, smiles2) in smartsref:
            patt1 = Chem.MolFromSmarts(smarts1)
            patt2 = Chem.MolFromSmiles(smiles2)
            if(m.HasSubstructMatch(patt1)):
                newmol = AllChem.ReplaceSubstructs(m, patt1, patt2)
                for ion in newmol:
                    i += 1
                    ion = Chem.MolToSmiles(ion,isomericSmiles=True)
                    if ion not in ionsmiles: ionsmiles.append(ion)

    number_sites = len(ionsmiles) - len(tautomers)
    for i in range(number_sites - 1):
        for ion in ionsmiles:
            m = Chem.MolFromSmiles(ion)
            for (smarts1, smiles2) in smartsref:
                patt1 = Chem.MolFromSmarts(smarts1)
                patt2 = Chem.MolFromSmiles(smiles2)
#               print ion, smarts1,smiles2
                if(m.HasSubstructMatch(patt1)):
                     newmol = AllChem.ReplaceSubstructs(m, patt1, patt2)
                     for new_ion in newmol:
                         smiles = Chem.MolToSmiles(new_ion,isomericSmiles=True)
#                        print "2",ion,smiles,smarts1,smiles2
                         if smiles not in ionsmiles: 
                             ionsmiles.append(smiles)
#                             print ion,smarts1,smiles2,smiles

#   ionsmiles.insert(0, neutral_smiles)

    for i,ion in enumerate(ionsmiles):
        charge = Chem.GetFormalCharge(Chem.MolFromSmiles(ion))
        if abs(charge) <= max_charge:
            m = Chem.MolFromSmiles(ion)
            charged_atoms = m.GetSubstructMatches(Chem.MolFromSmarts('[+,-]'))
            rings = m.GetRingInfo().AtomRings()
#           print ion,charged_atoms,in_diff_rings(m,charged_atoms,rings) 
            if len(charged_atoms) > 1 and len(rings) > 0:
               if in_diff_rings(m,charged_atoms,rings): 
                  print name+"_"+str(charge)+"="+str(i), ion
            else:
               print name+"_"+str(charge)+"="+str(i), ion

