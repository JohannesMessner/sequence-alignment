#####################################################
#######        COMPUTATIONAL BIOLOGY         ########
#######             HOMEWORK 1               ########
#####################################################
#                                                   #
# Implement the pairwise alignment algorithms       #
# Needleman-Wunsch and Smith-Waterman.              #
#                                                   #
#####################################################
#####################################################

# In all functions the following parameters are the same:
# seqA: the first sequence to align
# seqB: the second sequence to align
# score_gap: score for a gap
# score_match: score for a character match
# score_mismatch: score for a character mismatch
# local: True if alignment is local, False otherwise

init_score_matrix = function(nrow, ncol, local, score_gap) {
    # Initialize the score matrix with zeros.
    # If the alignment is global, the leftmost column and the top row will have incremental gap scores,
    # i.e. if the gap score is -2 and the number of columns is 4, the top row will be [0, -2, -4, -6].
    # nrow: number of rows in the matrix
    # ncol: number of columns in the matrix
    
    if(local) {
        score_matrix <- matrix(nrow = nrow, ncol = ncol, data = rep(0, nrow*ncol))
    } else {
        score_matrix <- matrix(nrow = nrow-1, ncol = ncol-1, data = rep(0, (nrow-1)*(ncol-1)))
        
        if (nrow <= 1){
            coldata <- score_gap
        } else {
            coldata <- c(1:(nrow-1)) * score_gap
        }
        
        col <- matrix(nrow = nrow-1, ncol=1, data = coldata)
        score_matrix <- cbind(col, score_matrix)
        
        rowdata <- c(0:(ncol-1)) * score_gap
        row <- matrix(nrow=1, ncol=ncol, data=rowdata)
        score_matrix <- rbind(row, score_matrix)
    }
  
    # Return the initialized empty score matrix
    # score_matrix: nrow by ncol matrix
    
    return(score_matrix)
}

init_path_matrix = function(nrow, ncol, local) {
    # Initialize the path matrix with empty values, except the top row and the leftmost column if global alignment.
    # If global alignment, the top row has "left" on all positions except 1st.
    # Similarly, leftmost column has "up" on all positions except 1st.
    # nrow: number of rows in the matrix
    # ncol: number of columns in the matrix

    if (local) {
        rowsymbol <- ""
        colsymbol <- ""
    } else {
        rowsymbol <- "up"
        colsymbol <- "left"
    }
    
    path_matrix <- matrix(nrow=nrow-1, ncol=ncol-1, data = rep("", (ncol-1)*(nrow-1)))
    
    col <- rep(rowsymbol, nrow-1)
    path_matrix <- cbind(col, path_matrix)
    
    row <- rep(colsymbol, ncol-1)
    row <- c("", row)
    path_matrix <- rbind(row, path_matrix)
    
    path_matrix <- unname(path_matrix)
    
    # Return the initialized empty path matrix
    # path_matrix: nrow by ncol matrix
    return(path_matrix)
}

get_best_score_and_path = function(row, col, nucA, nucB, score_matrix, score_gap, score_match, score_mismatch, local) {
    # Compute the score and the best path for a particular position in the score matrix
    # nucA: nucleotide in sequence A
    # nucB: nucleotide in sequence B
    # row: row-wise position in the matrix
    # col: column-wise position in the matrix
    # score_matrix: the score_matrix that is being filled out

    diagscore <- score_matrix[row-1, col-1]
    if (nucA == nucB) {
        match_mismatch_score = score_match
    } else {
        match_mismatch_score = score_mismatch
    }
    diagscore <- diagscore + match_mismatch_score
    leftscore <- score_matrix[row, col-1] + score_gap
    upscore <- score_matrix[row-1, col] + score_gap
    
    scores <- c(diagscore, leftscore, upscore)
    paths <- c("diag", "left", "up")
    names(paths) <- scores
    
    score <- max(scores)
    path <- paths[toString(score)]
    
    if(local && score < 0) {
        score <- 0
        path <- "-"
    }

    # Return the best score for the particular position in the score matrix
    # In the case that there are several equally good paths available, return any one of them.
    # score: best score at this position
    # path: path corresponding to the best score, one of ["diag", "up", "left"] in the global case and of ["diag", "up", "left", "-"] in the local case
    return(list("score"=score, "path"=path))
}

fill_matrices = function(seqA, seqB, score_gap, score_match, score_mismatch, local, score_matrix, path_matrix) {
    # Compute the full score and path matrices
    # score_matrix: initial matrix of the scores
    # path_matrix: initial matrix of paths
    
    if (nrow(score_matrix)>1 && ncol(score_matrix)>1){
        for (row in 2:nrow(score_matrix)){
            for (col in 2:ncol(score_matrix)){
                #seqB is assumed to be "on top" of the table, seqA to the left of it.
                nucA = substring(seqA, row-1, row-1)
                nucB = substring(seqB, col-1, col-1)
                score_and_path <- get_best_score_and_path(row=row, col=col, nucA=nucA, nucB=nucB,
                                                      score_matrix=score_matrix, score_gap = score_gap,
                                                      score_match = score_match, score_mismatch = score_mismatch,
                                                      local=local)

                score_matrix[row, col] <- score_and_path[["score"]]
                path_matrix[row, col] <- score_and_path[["path"]]
            }
        }
    }

    # Return the full score and path matrices
    # score_matrix: filled up matrix of the scores
    # path_matrix: filled up matrix of paths
    return(list("score_matrix"=score_matrix, "path_matrix"=path_matrix))
}

get_best_move = function(nucA, nucB, path, row, col) {
    # Compute the aligned characters at the given position in the score matrix and return the new position,
    # i.e. if the path is diagonal both the characters in seqA and seqB should be added,
    # if the path is up or left, there is a gap in one of the sequences.
    # nucA: nucleotide in sequence A
    # nucB: nucleotide in sequence B
    # path: best path pre-computed for the given position
    # row: row-wise position in the matrix
    # col: column-wise position in the matrix

    if (is.na(path)) {
        newrow <- NA
        newcol <- NA
        char1 <- nucA
        char2 <- nucB
    } else if (path == "diag") {
        newrow <- row - 1
        newcol <- col -1
        char1 <- nucA
        char2 <- nucB
    } else if (path == "up") {
        newrow <- row - 1
        newcol <- col
        char1 <- nucA
        char2 <- "-"
    } else if (path == "left") {
        newrow <- row
        newcol <- col - 1
        char1 <- "-"
        char2 <- nucB
    } else {
        newrow <- NA
        newcol <- NA
        char1 <- nucA
        char2 <- nucB
    }

    # Return the new row and column and the aligned characters
    # newrow: row if gap in seqA, row - 1 otherwise
    # newcol: col if gap in seqB, col - 1 otherwise
    # char1: '-' if gap in seqA, appropriate character if a match
    # char2: '-' if gap in seqB, appropriate character if a match
    return(list("newrow"=newrow, "newcol"=newcol, "char1"=char1, "char2"=char2))
}

get_best_alignment = function(seqA, seqB, score_matrix, path_matrix, local) {
    # Return the best alignment from the pre-computed score matrix
    # score_matrix: filled up matrix of the scores
    # path_matrix: filled up matrix of paths

    if(! local) {
        row <- nrow(path_matrix)
        col <- ncol(path_matrix)
    } else {
        max_points <- max(score_matrix)
        coordinates <- which(score_matrix == max_points, arr.ind = TRUE)[1,]
        row <- unname(coordinates["row"])
        col <- unname(coordinates["col"])
    }
    score <- score_matrix[row, col]
    
    alignment <- c("","")
    while (row > 0 && col > 0) {
        nucA <- substr(seqA, row-1, row-1)
        nucB <- substr(seqB, col-1, col-1)
        path <- path_matrix[row, col]
        
        move <- get_best_move(nucA = nucA, nucB = nucB, path = path, row=row, col=col)
        
        row <- move[["newrow"]]
        col <- move[["newcol"]]
        char1 <- move["char1"]
        char2 <- move["char2"]
        
        if (is.na(row) || is.na(col)) {
            break
        } else {
            alignment[1] <- paste0(char1, alignment[1])
            alignment[2] <- paste0(char2, alignment[2])
        }
    }
    # Return the best score and alignment (or one thereof if there are multiple with equal score)
    # score: score of the best alignment
    # alignment: the actual alignment in the form of a vector of two strings
    return(list("score"=score, "alignment"=alignment))
}

align = function(seqA, seqB, score_gap, score_match, score_mismatch, local) {
    # Align the two sequences given the scoring scheme
    # For testing purposes, use seqA for the rows and seqB for the columns of the matrices
  
    # Initialize score and path matrices
    nrow <- nchar(seqA) + 1
    ncol <- nchar(seqB) + 1
    score_matrix <- init_score_matrix(nrow = nrow, ncol = ncol, local = local, score_gap = score_gap)
    
    path_matrix <- init_path_matrix(nrow = nrow, ncol = ncol, local = local)
    
    # Fill in the matrices with scores and paths using dynamic programming
    
    matrices <- fill_matrices(seqA = seqA, seqB = seqB, score_gap = score_gap, score_match = score_match,
                              score_mismatch = score_mismatch, local = local, score_matrix = score_matrix,
                              path_matrix = path_matrix)
    score_matrix <- matrices[["score_matrix"]]
    
    path_matrix <- matrices[["path_matrix"]]
  
    # Get the best score and alignment (or one thereof if there are multiple with equal score)
    
    result <- get_best_alignment(seqA = seqA, seqB = seqB, score_matrix = score_matrix,
                                 path_matrix = path_matrix, local = local)
    
    # Return the best score and alignment (or one thereof if there are multiple with equal score)
    # Returns the same value types as get_best_alignment
    return(result)
}

test_align = function() {
    seqA = "TCACACTAC"
    seqB = "AGCACAC"
    score_gap = -2
    score_match = +3
    score_mismatch = -1
    local = F
    result = align(seqA, seqB, score_gap, score_match, score_mismatch, local)
    print(result$alignment)
    print(result$score)
}

test_align()