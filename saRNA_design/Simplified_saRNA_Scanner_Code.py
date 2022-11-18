from pylab import *

promotername = input('Promotername:')

AAE = -4.26
ATE = -3.67
TAE = -2.5
CAE = -6.12
GTE = -6.09
CTE = -5.4
GAE = -5.51
CGE = -9.07
GCE = -9.36
GGE = -7.66

fudge = 0.00000001

end = 4

stab = 0.01

ind18 = 1000.

ind7 = 1000.

flank = 1.

filepath = 'saRNAs_design'


def checksequence(sequence):
    sequence = sequence.upper()

    i = 0
    flag = True
    while i < len(sequence) and flag:
        if sequence[i] != 'A' and sequence[i] != 'G' and sequence[i] != 'C' and sequence[i] != 'T' and sequence[
            i] != 'U':
            flag = False
        i = i + 1

    if flag:
        sequence = sequence.replace('T', 'U')
    else:
        sequence = 'failure'

    return sequence


def maketargets(sequence):
    i = 0
    targets = []
    while i < len(sequence) - 18:
        targets = targets + [sequence[i:i + 19]]
        i = i + 1
    return targets


def countGC(target):
    i = 0
    GC = 0
    while i < len(target):
        if target[i] == 'G' or target[i] == 'C':
            GC = GC + 1
        i = i + 1

    GC = float(GC)
    GC = GC / len(target)
    return GC



def checkconsecutive(sequence, position):
    i = 4
    flag = False
    while i < 19 and flag == False:
        if sequence[position - 4 + i] == sequence[position - 3 + i] and sequence[position - 3 + i] == sequence[
            position - 2 + i] and sequence[position - 2 + i] == sequence[position - 1 + i] and sequence[
            position - 1 + i] == sequence[position + i]:
            flag = True

        i = i + 1

    return flag



def endbase(target):
    if target[0] == 'A' or target[0] == 'U':
        energy = AAE / 2
    elif target[0] == 'G' or target[0] == 'C':
        energy = GGE / 2
    else:
        print("")
        print("Error: none of the standard bases were detected.")

    E5prime = energy

    if target[-1] == 'A' or target[-1] == 'U':
        energy = AAE / 2
    elif target[-1] == 'G' or target[-1] == 'C':
        energy = GGE / 2
    else:
        print("")
        print("Error: none of the stardard bases were detected.")

    E3prime = energy

    if E3prime - E5prime > fudge:
        recog = True
    else:
        recog = False

    return recog



def endstability(target):
    i = 0
    E5prime = 0.
    while i < end - 1:
        if target[i:i + 2] == 'AA' or target[i:i + 2] == 'UU':
            energy = AAE
        elif target[i:i + 2] == 'AU':
            energy = ATE
        elif target[i:i + 2] == 'UA':
            energy = TAE
        elif target[i:i + 2] == 'CA' or target[i:i + 2] == 'AC':
            energy = CAE
        elif target[i:i + 2] == 'GU' or target[i:i + 2] == 'UG':
            energy = GTE
        elif target[i:i + 2] == 'CU' or target[i:i + 2] == 'UC':
            energy = CTE
        elif target[i:i + 2] == 'GA' or target[i:i + 2] == 'AG':
            energy = GAE
        elif target[i:i + 2] == 'CG':
            energy = CGE
        elif target[i:i + 2] == 'GC':
            energy = GCE
        elif target[i:i + 2] == 'GG' or target[i:i + 2] == 'CC':
            energy = GGE
        else:
            print("")
            print("Error: none of the recognized nearest-neighbor pairings were detected.")

        E5prime = E5prime + energy
        i = i + 1

    i = len(target) - end
    E3prime = 0.
    while i < len(target) - 1:
        if target[i:i + 2] == 'AA' or target[i:i + 2] == 'UU':
            energy = AAE
        elif target[i:i + 2] == 'AU':
            energy = ATE
        elif target[i:i + 2] == 'UA':
            energy = TAE
        elif target[i:i + 2] == 'CA' or target[i:i + 2] == 'AC':
            energy = CAE
        elif target[i:i + 2] == 'GU' or target[i:i + 2] == 'UG':
            energy = GTE
        elif target[i:i + 2] == 'CU' or target[i:i + 2] == 'UC':
            energy = CTE
        elif target[i:i + 2] == 'GA' or target[i:i + 2] == 'AG':
            energy = GAE
        elif target[i:i + 2] == 'CG':
            energy = CGE
        elif target[i:i + 2] == 'GC':
            energy = GCE
        elif target[i:i + 2] == 'GG' or target[i:i + 2] == 'CC':
            energy = GGE
        else:
            print("")
            print("Error: none of the recognized nearest-neighbor pairings were detected.")

        E3prime = E3prime + energy
        i = i + 1

    stabilitydiff = E3prime - E5prime
    return stabilitydiff



def checkpos19(target):
    flag = False
    if target[18] == 'A':
        flag = True
    return flag



def checkpos18(target):
    indicator = 0
    if target[17] == 'A':
        indicator = 2
    elif target[17] == 'U':
        indicator = 1
    return indicator



def checkpos7(target):
    indicator = 0
    if target[6] == 'U':
        indicator = 1
    return indicator



def checkflanks(sequence, position):
    i = position + 19
    score = 0
    while i < position + 23:
        if sequence[i] == 'A' or sequence[i] == 'U':
            score = score + 1
        i = i + 1
    return score



def CpGislanddetector(sequence):
    if len(sequence) > 200:

        CpGverdicts = []
        GC = (sequence.count('G', 0, 200) + sequence.count('C', 0, 200)) / 200.
        CpG = sequence.count('CG', 0, 200)

        CpGexp = (sequence.count('C', 0, 200) * sequence.count('G', 0, 200)) / 200.

        if GC > 0.5 and CpG > 0.6 * CpGexp and CpG > 7:
            CpGverdicts = CpGverdicts + [True]
        else:
            CpGverdicts = CpGverdicts + [False]

        i = 1

        while i < len(sequence) - 199:
            GCchange = 0
            if sequence[i - 1] == 'G' or sequence[i - 1] == 'C':
                GCchange = GCchange - 1
            if sequence[i + 199] == 'G' or sequence[i + 199] == 'C':
                GCchange = GCchange + 1
            GC = ((GC * 200.) + GCchange) / 200.
            CpGexp = (sequence.count('C', i, i + 200) * sequence.count('G', i, i + 200)) / 200.

            CpGchange = 0
            if sequence[i - 1] == 'C' and sequence[i] == 'G':
                CpGchange = CpGchange - 1
            if sequence[i + 198] == 'C' and sequence[i + 199] == 'G':
                CpGchange = CpGchange + 1
            CpG = CpG + CpGchange

            if GC > 0.5 and CpG > 0.6 * CpGexp and CpG > 7:
                CpGverdicts = CpGverdicts + [True]
            else:
                CpGverdicts = CpGverdicts + [False]

            i = i + 1

        rawislandboundaries = []
        i = 0
        j = 0
        while i < len(CpGverdicts):
            if CpGverdicts[i]:
                rawislandboundaries = rawislandboundaries + [[]]
                rawislandboundaries[j] = rawislandboundaries[j] + [i, i + 199]
                j = j + 1
            i = i + 1

        if rawislandboundaries == []:
            islandboundaries = []
            print('No CpG island regions were identified.')
        else:
            islandboundaries = [rawislandboundaries[0]]
            i = 0
            j = 0
            while i < len(rawislandboundaries):
                if rawislandboundaries[i][0] < islandboundaries[j][1]:
                    islandboundaries[j][1] = rawislandboundaries[i][1]

                else:
                    j = j + 1
                    islandboundaries = islandboundaries + [rawislandboundaries[i]]
                i = i + 1
    else:
        print("")
        print(
            "Warning: the submitted sequence is less than 200 bp long. This makes positive CpG island identification impossible.")

    return islandboundaries



def inisland(targetposition, islandboundaries):
    i = 0
    flag = False
    while i < len(islandboundaries) and flag == False:
        if (targetposition >= islandboundaries[i][0] and targetposition < islandboundaries[i][1]) or (
                targetposition + 18 > islandboundaries[i][0] and targetposition + 18 <= islandboundaries[i][1]):
            flag = True
        i = i + 1

    return flag


def CpGsitedetector(sequence):
    i = 1
    CpGs = []
    while i < len(sequence):
        if sequence[i - 1:i + 1] == 'CG':
            CpGs = CpGs + [i - 1]
        i = i + 1

    return CpGs


def insite(targetposition, CpGs):
    i = 0
    flag = False
    while i < len(CpGs) and flag == False:
        if (CpGs[i] >= targetposition and CpGs[i] <= targetposition + 18) or (
                CpGs[i] + 1 >= targetposition and CpGs[i] + 1 <= targetposition + 18):
            flag = True
        i = i + 1

    return flag


def complement(target):
    i = 0
    comp = ''
    while i < len(target):
        if target[i] == 'A':
            comp = comp + 'U'
        elif target[i] == 'G':
            comp = comp + 'C'
        elif target[i] == 'C':
            comp = comp + 'G'
        elif target[i] == 'U':
            comp = comp + 'A'
        else:
            print("")
            print("An error occurred taking the reverse complement of the sequence.")

        i = i + 1

    return comp


sequence = input('Input promoter sequence of interest:')
print("")
print("Checking sequence..")
sequence = checksequence(sequence)

print("Scanning for CpG elements..")
CpGislands = CpGislanddetector(sequence)
CpGsites = CpGsitedetector(sequence)

if sequence == 'failure':
    print("")
    print(
        "Error: your sequence is in an invalid format. Only the characters 'A', 'G', 'C', and 'T' or 'U' (uppercase or lowercase)) are allowed.")
else:
    print("Generating list of potential targets..")
    alltargets = maketargets(sequence)
    i = 3
    survivors = []
    print("Validating targets..")
    while i < len(alltargets) - 4:
        flag = True

        if checkconsecutive(sequence, i):
            flag = False

        if endstability(alltargets[i]) < 0:
            flag = False
        # default min 0.4 max 0.65
        if countGC(alltargets[i]) < 0.4 or countGC(alltargets[i]) > 0.65:
            flag = False

        if checkpos19(alltargets[i]) == False:
            flag = False

        if checkpos18(alltargets[i]) == 0:
            flag = False

        if endbase(alltargets[i]) == False:
            flag = False

        # if inisland(i, CpGislands):
        #    flag = False

        if insite(i, CpGsites):
            flag = False

        if flag:
            survivors = survivors + [[i, alltargets[i]]]

        i = i + 1

    i = 0
    survivorscores = []
    print("Ranking target sequences..")
    while i < len(survivors):
        stabilitydiff = endstability(survivors[i][1])
        indicator18 = checkpos18(survivors[i][1])
        indicator7 = checkpos7(survivors[i][1])
        flankscore = checkflanks(sequence, survivors[i][0])
        score = stab * stabilitydiff + ind18 * indicator18 + ind7 * indicator7 + flank * flankscore
        survivorscores = survivorscores + [score]

        i = i + 1

    if len(survivors) != len(survivorscores):
        print("Error: scoring proceeded incorrectly for an unknown reason.")

    i = 0
    while i < len(survivors):
        comp = complement(survivors[i][1])
        survivors[i] = survivors[i] + [survivorscores[i]] + [comp]
        i = i + 1

    sortedsurvivors = sorted(survivors, key=lambda survivors: survivors[2])

    print("Writing output to file..")
    output = open(filepath + 'dsRNA_Candidates_Promoter_' + promotername + '.txt', 'w')

    if len(sortedsurvivors) == 0:
        output.write('No valid targets were detected in this promoter.')

    i = len(sortedsurvivors) - 1
    while i > -1:
        output.write('-------------------- dsRNA ' + str(len(sortedsurvivors) - i) + ' ---------------------\n\n')
        output.write('Position on target strand: ' + str(sortedsurvivors[i][0] + 1) + '\n')
        output.write('Overall score: ' + str(sortedsurvivors[i][2]) + '\n\n')
        output.write('Guide strand         3\' UU' + sortedsurvivors[i][3] + '    5\'\n')
        output.write('                          |||||||||||||||||||\n')
        output.write('Passenger strand     5\'   ' + sortedsurvivors[i][1] + 'UU  3\'\n\n\n')
        i = i - 1

    output.close()

print("")
print(
    "dsRNAs have been generated and ranked according to scores for optimality. Sequences are available in file dsRNA_Candidates_Promoter_" + promotername + ".txt.")
print("End of program.")

