import sys

#input of the two reads
first = input('Insert the first sequence: ')
second = input('Insert the second sequence: ')

#defining the alphabet for the substitution matrix
alphabet = "ATCG"

#running: python GiuliaGrotto.py scoreMatch scoreMismatch scoreGap  
#match = int(sys.argv[1])
#mismatch = int(sys.argv[2])
#gapPenalty = int(sys.argv[3])

#score assigned to each possible situation
match = 3
mismatch = -3
gapPenalty = -2
gapAffine = 1

#-MAIN

def main(): 

    #calculating score matrix dimensions'
    rows = len(first) + 1
    cols = len(second) + 1

    #fill in substitution matrix
    subMatrix = buildSubMatrix(alphabet, match, mismatch)

    #creating and filling in the score Matrix and the path Matrix 
    scoreMatrix, pathMatrix, maxScore, maxPos = createScoreMatrix(subMatrix, rows, cols)

    print("Score matrix: \n")
    print_matrix(scoreMatrix)
    print("Max score: ", maxScore)
    print("In position: ", maxPos)

    #list with the position of the best scores (>50% wrt maximum score) with decreasing order (HIGHER TO LOWER as requested) in score matrix
    bestScores = findBest(scoreMatrix, rows, cols, maxScore, 70) 

    #print all the alignment starting from the best scores calculated above, without overlaps
    for i, j, k in bestScores:
        #visited = pathMatrix[i][j][2]

        #if the start position of the alignment is not inside another path better than this
        #if not visited:
        #reconstructing the best alignment of the two reads
        alignment, nmatches, ngaps, nmismatches = traceback(pathMatrix, scoreMatrix, [i, j])

        #this condition check if there is only one gap and if there are at least two occurrences of at least three mismatches
        if ngaps == 1 and condition(alignment[2][::-1]):
            print("\nThe alignment is: \n" + alignment[1][::-1] + "\n" + alignment[2][::-1] + "\n" + alignment[0][::-1])
            print("The score of the alignment is: ", k)
            print("The alignment length is: ", len(alignment[0]))
            print('With {0} matches, {1} mismatches and {2} gaps'.format(nmatches, nmismatches, ngaps))

    return(0)  


#--FUNCTIONS  

#create substitution matrix for the given alphabet
def buildSubMatrix(alphabet, match, mismatch):
    subMatrix = [[0 for x in range(len(alphabet))] for y in range(len(alphabet))]
    
    for i in range(len(alphabet)):
            for j in range(len(alphabet)):
                if alphabet[i] == alphabet[j] :
                    subMatrix[i][j] = match
                elif alphabet[i] != alphabet[j]:
                    subMatrix[i][j] = mismatch

    return subMatrix


#calculate the score if match or mismatch from substitution matrix
def subScore(subMatrix, letter1, letter2):
    a = alphabet.index(letter1)
    b = alphabet.index(letter2)

    return subMatrix[a][b]


#calculate the penalty for the gap (affine or linear)
def gap_penalty(k) : 
    #return (gapPenalty*k) + gapAffine #affine
    return gapPenalty*k #linear
    

#find the index of the maximum score in the list containing scores for n gaps
def maxGap(gapValues):
    maxGap = max(gapValues)
    return gapValues.index(maxGap)


#create and fill in the score matrix and the path matrix used for the traceback    
def createScoreMatrix(subMatrix, rows, cols):
    scoreMatrix = [[0 for x in range(cols)] for y in range(rows)] 
    #a cell of the path matrix will contain:
    #-the number of steps, 
    #-the direction where the score come from,
    #-if the node has been visited or not during the traceback
    pathMatrix = [[[0, 'NULL', False] for x in range(cols)] for y in range(rows)]

    maxScore = 0
    maxPos = (0, 0)

    #fill in score matrix
    for x in range(1, rows):
        for y in range(1, cols):

            #calculating the score for every possible number of gaps til the current position (from up and left)
            same_row = [(scoreMatrix[x][y-l]+gap_penalty(l)) for l in range(1,y+1)]
            same_col = [(scoreMatrix[x-k][y]+gap_penalty(k)) for k in range(1,x+1)]
    
            #take the maximum score between them 
            upscore = max(same_col)
            leftscore = max(same_row)

            #score from the substitution matrix if match or mismatch
            alike = subScore(subMatrix, first[x-1], second[y-1])
            diagscore = scoreMatrix[x-1][y-1] + alike

            #take the index, as the number of gaps corresponding to the best score from up or left 
            maxGapUp = maxGap(same_col)
            maxGapLeft = maxGap(same_row)

            score = max(0, upscore, diagscore, leftscore)
            scoreMatrix[x][y] = score

            if score >= maxScore:
                maxScore = score
                maxPos = (x, y)

            #saving ONE direction where the score is coming from with the number of steps
            #U: up
            #L: left
            #D: diagonal
            direction = 'NULL'
            if score == upscore:
                direction = [maxGapUp + 1, 'U', False]
            elif score == leftscore:
                direction = [maxGapLeft + 1, 'L', False]
            elif score == diagscore:
                if alike > 0:
                    direction = [1, 'D', False]
                else:
                    direction = [-1, 'D', False]

            pathMatrix[x][y] = direction

    return scoreMatrix, pathMatrix, maxScore, maxPos

#find the score greater than the percentage of the maximum score and their position, with score decreasing order
def findBest(scoreMatrix, rows, cols, maxScore, percentage):
    least = maxScore*percentage/100
    bestScores = []

    for i in range(1, rows):
        for j in range(1, cols):
            if scoreMatrix[i][j] > least :
                bestScores.append([i, j, scoreMatrix[i][j]])

    bestScores = sorted(bestScores, key = lambda x: x[2])
    return bestScores[::-1]

#traceback which reconstruct the best local alignment
def traceback(pathMatrix, scoreMatrix, maxPos):
    i = maxPos[0]
    j = maxPos[1]

    #save the two strings for the alignment and a graphic
    alignment = ['', '', '']

    ngaps, nmatches, nmismatches = 0, 0, 0

    while scoreMatrix[i][j] != 0:
        n = pathMatrix[i][j][0]
        direction = pathMatrix[i][j][1]
        #set visited as true
        pathMatrix[i][j][2] = True
        #print(direction)
        #print("Score: ", scoreMatrix[i][j])
        if direction == 'D':
            alignment[0] = alignment[0] + first[i-1]
            alignment[1] = alignment[1] + second[j-1]
            pathMatrix[i][j][2] = True
            if n == 1:
                nmatches += 1
                alignment[2] = alignment[2] + '|'
            else:
                nmismatches += 1
                alignment[2] = alignment[2] + ':'
            i = i - 1
            j = j - 1
        elif direction == 'U':
            for k in range(n):
                alignment[0] = alignment[0] + '-'
                alignment[1] = alignment[1] + first[i-1]
                alignment[2] = alignment[2] + ' '
                pathMatrix[i][j][2] = True
                ngaps += 1
                i = i - 1
        elif direction == 'L':
            for l in range(n):
                alignment[0] = alignment[0] + second[j-1]
                alignment[1] = alignment[1] + '-'
                alignment[2] = alignment[2] + ' '
                pathMatrix[i][j][2] = True
                ngaps += 1
                j = j - 1
        
    return alignment, nmatches, ngaps, nmismatches

#function that find if there are 2 occurrences of 3 consecutive matches
#using my past code I decided to use the graphics that I made for the alignment output, so that I don't have to intefer with the traceback code
#In the string for every step there is "|" for a match, ":" for a mismatch and " " for a gap
#So, if there are 2 occurrences of "|||", means that there are 2 occurrences of 3 consecutive matches
#I could use another string with like a character for every type, for example M: matches, m:mismathces and g:gap, but it was worthless to save another string 
#having another one with the same meaning
def condition(stringa):
    isRight = True
    for i in range(2):
        index = stringa.find("|||")
        if index != -1:
            stringa = stringa[(index + 3):len(stringa)]
            isRight = isRight and True
        else:
            isRight = isRight and False
    return isRight
    
# graphical display of matrix 
def print_matrix(matrix): 
    print('\n'.join([''.join(['     {:4}'.format(item) for item in row])for row in matrix]))


if __name__ == '__main__':
    sys.exit(main())
    