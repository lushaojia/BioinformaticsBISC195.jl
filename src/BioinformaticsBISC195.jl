module BioinformaticsBISC195

export normalizeDNA,
       composition,
       gc_content,
       complement,
       reverse_complement,
       parse_fasta,
       removeUnder25k,
       uniqueKmers,
       compareKmerSets,
       nwscore,
       nwsetupmatrix,
       nwscorematrix,
       nwalign,
       alignmentHeatmap,
       mis_matchSeq,
       AlignmentSubsequences,
       mis_matchSubsequences302
       


# # uncomment the following line if you intend to use BioSequences types
# using BioSequences

"""
    normalizeDNA(::AbstractString)

Ensures that a sequence only contains valid bases
(or `'N'` for unknown bases). 

Ambiguous bases other
than 'N' (i.e. one of "RYSWKMBDHV.-") are changed to 'N'.

Returns a String of a DNA sequence with valid bases only.
Bases are all uppercase and will be one of `'A', 'G', 'C', 'T', or 'N' (ambiguous base)`
"""
function normalizeDNA(seq)
    seq = uppercase(string(seq))
    seqStr = ""
    
    for base in seq
        if occursin(base, "RYSWKMBDHV.-")
            seqStr *= 'N'
        else
            seqStr *= base
        end
        # note: `N` indicates an unknown base
        occursin(seqStr[length(seqStr)], "AGCTN") || error("Invalid base, $base")
    end
    return seqStr # change to `return LongDNASeq(seq)` if you want to try to use BioSequences types
end

"""
    composition(sequence)

Counts the number of each type of base
in a DNA sequence and returns a tuple with the counts in this order:
1. # of 'A's  
2. # of 'G's 
3. # of 'C's 
4. # of 'T's
5. # of 'N's

Throws an error when an invalid base (not 'A', 'C',
'G', 'T', or 'N') is encountered.

Examples  
≡≡≡≡≡≡≡≡≡≡
    julia> composition("ACCGGGTTTTN")
        
    

    julia> composition("AAX")
    ERROR: Invalid base, X
"""
function composition(sequence)
    #normalizedSeq = normalizeDNA(sequence)
    return count("A", sequence), count("G", sequence), count("C", sequence), count("T", sequence), count("N", sequence)
end

"""
    gc_content(seq)

Calculates the GC ratio (the total number of G and C bases 
divided by the total length of the sequence) of a DNA sequence.

Accomodates sequences with ambiguous bases (e.g. bases R, Y, S, W, K, M, B, D, H, V, N and gaps . or -).

Examples  
≡≡≡≡≡≡≡≡≡≡

    julia> gc_content("AATG")
    0.25

    julia> gc_content("cccggg") * 100
    100.0

    julia> gc_content("ATta")
    0.0

    julia> gc_content("ccccggggn")
    0.8888888888888888

    julia> gc_content("ATty")
    Error: Invalid base Y encountered
"""
function gc_content(seq)
    baseCount = composition(seq)
    return (baseCount[2] + baseCount[3])/(length(seq))
end

"""
    complement(sequence)

Returns the complementary strand of `sequence`.
The complementary base pairings are as follows:

    A <-> T
    G <-> C
    N <-> N

Accepts uppercase or lowercase `String`,
but always returns an uppercase `String`.

Accomodates sequences with ambiguous bases (e.g. bases R, Y, S, W, K, M, B, D, H, V, N and gaps . or -).

If a valid base is not provided, the function throws an error.
"""

function complement(sequence)
    compDict = Dict('A'=>'T', 'T'=>'A', 'G'=>'C', 'C'=>'G', 'N'=>'N')
    normalizedSeq = normalizeDNA(sequence)

    compStr = ""
    for base in normalizedSeq
        !(base in keys(compDict)) && error("Invalid base, $base")
        compStr *= compDict[base]
    end
    return compStr
end

"""
    reverse_complement(sequence)

Takes a DNA sequence and returns the reverse complement
of that sequence.

Takes lowercase or uppercase sequences,
but always returns uppercase.

Accomodates sequences with ambiguous bases (e.g. bases R, Y, S, W, K, M, B, D, H, V, N and gaps . or -).

Examples
≡≡≡≡≡≡≡≡≡≡
    julia> reverse_complement("AAATTT")
    "AAATTT"

    julia> reverse_complement("GCAT")
    "ATGC"

    julia> rc = reverse_complement("TTGGG");

    julia> println(rc)
    CCCAA

    julia> reverse_complement("ATTAGC")
    "GCTAAT"

    julia> reverse_complement("ATN")
    "NAT"
"""

function reverse_complement(sequence)
    complSeq = complement(sequence)
    revComplSeq = reverse(complSeq)
    return revComplSeq
end


"""
    function parse_fasta(path; DNA=true)

Reads a fasta-formated file and returns 2 vectors,
one containing parsed headers 
(which do not contain the character '>' at the beginning),
the other containing the entire sequence as a `String`.

If the kwarg DNA is true, the function validates sequences to see if they contain DNA bases only. Accomodates 
sequences with ambiguous bases (e.g. bases R, Y, S, W, K, M, B, D, H, V, N and gaps . or -).

If the kwarg DNA is false, the function does not validate the sequences to see if they contain DNA bases only.

Examples
≡≡≡≡≡≡≡≡≡
    julia> ex1 = parse_fasta("data/ex1.fasta");

    julia> ex1[1]
    2-element Array{String,1}:
    "ex1.1 | easy"
    "ex1.2 | multiline"

    julia> ex1[2]
    2-element Array{String,1}:
    "AATTATAGC"
    "CGCCCCCCAGTCGGATT"

    julia> ex2 = parse_fasta("data/ex2.fasta");
    ERROR: invalid base H
"""

function parse_fasta(path; DNA=true)
    headerArr = []
    seqStr = ""
    seqArr = []
    counter = 1

    for line in readlines(path)
        if startswith(line, '>')
            alteredHeader = line[2:length(line)]
            push!(headerArr, alteredHeader)
            if counter > 1
                # every time a header is encountered
                # pushes the concatenated DNA sequences that came before the header to seqArr
                # and
                # seqStr is re-set to an empty String
                push!(seqArr, seqStr)
                seqStr = ""
            end
            counter += 1
        else
            # for base in line
            #     occursin(base, "AGCTRYSWKMBDHVN.-") || error("Invalid base, $base")
            # end
            if DNA
                normalLine = normalizeDNA(line)
                seqStr *= normalLine
            else
                seqStr *= line
            end
        end
    end

    # push the last conactenated DNA sequence to seqArr
    push!(seqArr, seqStr)

    return headerArr, seqArr
end

"""
    function removeUnder25k(parsed_data)

Takes an argument with same data structure as the return value of parse_fasta(path).
(i.e. Takes an Array with 2 nested Arrays. The first nested Array contains `Strings` of
parsed headers roughly in the following format: "NC_045512.2 | China| Homo sapiens". 
The second nested Array contains `Strings` of entire DNA sequences whose base composition is strictly
AGCTN.)

Returns a value in the same data structure as described above, except that all DNA sequences shorter
than 25k basepairs and their headers are taken out.

Examples
≡≡≡≡≡≡≡≡≡
    julia> parsed_ex1 = parse_fasta("ex1.fasta")
    julia> removeUnder25k(parsed_ex1)
    2-element Vector{Vector{Any}}:
    []
    []

    julia> cleaned_data = removeUnder25k("/mnt/c/Users/lusha/Documents/Wellesley/BISC195/CoV2-genome-analysis/SARS-related-coronaviruses-genomes.fasta")
    julia> findall(seq -> length(seq)<25000, cleaned_data[2])
    Int64[]
"""
function removeUnder25k(parsed_data)
    deleteIndArr = findall(seq -> length(seq)<25000, parsed_data[2])
    deleteat!(parsed_data[1], deleteIndArr) # deleting the headers all at once
    deleteat!(parsed_data[2], deleteIndArr) # deleting the sequences all at once
    return parsed_data
end

"""
    function uniqueKmers(seq, k)

Takes a DNA sequence and returns an Array of all of its 
unique kmers of length k with no ambiguous bases. 
(Assumes that this DNA sequence has been processed with normalizeDNA(seq)
such that all of its non-N ambiguous bases have been converted to N's or that
this DNA sequence only contains N's as ambiguous bases.)  

Throws an error if k is longer than the length of the DNA sequence.

Examples
≡≡≡≡≡≡≡≡≡
    julia> uniqueKmers("ATATATGC", 2)
    4-element Vector{Any}: "AT" "TA" "TG" "GC"

    julia> uniqueKmers("GTNT", 2)
    1-element Vector{Any}: "GT"

    julia> uniqueKmers("AGTTTTATN", 17)
    k must be a positive integer less than the length of the sequence
"""
function uniqueKmers(seq, k)
    1 <= k <= length(seq) || error("k must be a positive integer less than the length of the sequence")
    stopindex = length(seq)-k+1
    uniqueKmersArr = []
    
    for i in 1:stopindex
        kmer = seq[i:i+k-1]
        (!(kmer in uniqueKmersArr) && !(occursin("N", kmer))) && push!(uniqueKmersArr, kmer)
    end 

    return uniqueKmersArr
end

"""
    function compareKmerSets(kmerSet1, kmerSet2)

Takes two unique kmer sets as Arrays. Returns their distance metric as a Float.
0<= distance metric <= 1. Two completely different kmer sets generate a distance metric of 1, while two
completely identical kmer sets generate a distance metric of 0. 
The distance metric is measured as one minus the number of kmers in the intersection of the two kmer sets
divided by the number of kmers in the union of the two kmer sets.

Examples
≡≡≡≡≡≡≡≡≡
    julia> kset1 = uniqueKmers("ATATG", 2)
    julia> kset2 = uniqueKmers("ATATG", 2)
    julia> compareKmerSets(kset1, kset2)
    0.0

    julia> typeof(compareKmerSets(kset1, kset2))
    Float64

    julia> compareKmerSets(["ATG", "CTA", "TAT"], ["ATG", "CTT", "TAT"])
    0.5 
"""

function compareKmerSets(kmerSet1, kmerSet2)
    # returns the distance metric
    return 1 - (length(intersect(kmerSet1, kmerSet2))/length(union(kmerSet1, kmerSet2)))
end

"""
    nwscore(base1, base2; mismatch = -1, match = 1)

Compares nucleotides for genomic sequences or
amino acids for protein sequences to see if they are 
matching=1 or mismatching=-1.

Examples
≡≡≡≡≡≡≡≡≡
    julia> nwscore('A','G')
    -1
    
    julia> nwscore('V','V')
    1

    julia> nwscore(nothing,nothing)
    ERROR: ArgumentError: Score for two gaps is not defined
    Stacktrace:
    [1] nwscore(#unused#::Nothing, #unused#::Nothing)
    @ Main ~/Documents/my_repo/BISC195Labs/src/needleman_wunch.jl:25
    [2] top-level scope
    @ REPL[81]:1
"""
function nwscore(base1::Char, base2::Char; match = 1, mismatch = -1)
    if base1 == base2
        return match
    else base1 != base2
        return mismatch
    end
end

function nwscore(::Nothing, ::Nothing)
    throw(ArgumentError("Score for two gaps is not defined"))
end

"""
    function nwsetupmatrix(s1, s2; gap=-1)

Takes two biological sequences and a gap score as a kwarg. 
Does not check the sequences for residue validity.
Returns a matrix where 
matrix[1, 1] is a Floating zero,
2 <= i <= length(s1) and 2 <= j <= length(s2) are Floating zeros, 
and the first row/column are headers of incrementing gap scores going right/down.

Examples
≡≡≡≡≡≡≡≡≡
    julia> nwsetupmatrix("MFV", "NVVIKV"; gap = -4)
    4×7 Matrix{Float64}:
    0.0  -4.0  -8.0  -12.0  -16.0  -20.0  -24.0
   -4.0   0.0   0.0    0.0    0.0    0.0    0.0
   -8.0   0.0   0.0    0.0    0.0    0.0    0.0
  -12.0   0.0   0.0    0.0    0.0    0.0    0.0
"""
function nwsetupmatrix(s1, s2; gap=-1)
    #setting up a matrix of all zero's
    setupmatrix = zeros(length(s1)+1, length(s2)+1)

    #making the horizontal and vertical "headers" with the appropriate gap scores
    for j in 2:length(s2)+1
        setupmatrix[1, j] = (j-1)*gap
    end
    for i in 2:length(s1)+1
        setupmatrix[i, 1] = (i-1)*gap
    end
    return setupmatrix
end

"""
    function nwscorematrix(s1, s2; match=1, mismatch=-1, gap=-1)

Takes two biological sequences as position arguments, and a match, mismatch, and gap score as kwargs.
Returns a dictionary with s1*s2 number of key=>value pairs, where keys are indices i, j of a matrix cell
as Tuples and values are Arrays containing the score in the matrix cell and the direction from which it came from (one of "above", "left", or "diagonal").

Examples
≡≡≡≡≡≡≡≡≡
    julia> nwscorematrix("AAATG", "AATG"; gap=-2, match=2)
    Dict{Any, Any} with 20 entries:
    (4, 5) => Any[1.0, "left"]
    (2, 5) => Any[-4.0, "left"]
    (6, 2) => Any[-6.0, "above"]
    (6, 3) => Any[-2.0, "above"]
    (6, 4) => Any[2.0, "above"]
    (5, 5) => Any[2.0, "left"]
    (3, 2) => Any[0.0, "above"]
    (3, 3) => Any[4.0, "diagonal"]
    (3, 4) => Any[2.0, "left"]
    (6, 5) => Any[6.0, "diagonal"]
    (4, 2) => Any[-2.0, "above"]
    (2, 2) => Any[2.0, "diagonal"]
    (4, 3) => Any[2.0, "above"]
    (2, 3) => Any[0.0, "left"]
    (3, 5) => Any[0.0, "left"]
    (4, 4) => Any[3.0, "diagonal"]
    (2, 4) => Any[-2.0, "left"]
    (5, 2) => Any[-4.0, "above"]
    (5, 3) => Any[0.0, "above"]
    (5, 4) => Any[4.0, "diagonal"]
"""
function nwscorematrix(s1, s2; match=1, mismatch=-1, gap=-1)
    scoremat = nwsetupmatrix(s1, s2; gap=gap)
    scoreDict = Dict()

    for i in 2:size(scoremat, 1)
        for j in 2:size(scoremat, 2)
            above = gap + scoremat[i-1, j]
            left = gap + scoremat[i, j-1]
            diagonal = nwscore(s1[i-1], s2[j-1]; match=match, mismatch=mismatch) + scoremat[i-1, j-1]
            scoreArr = [above, left, diagonal]
            maxScore = maximum(scoreArr)
            scoremat[i, j] = maxScore

            #pushing to scoreDict: keys are indices of current matrix cell in a Tuple; values are Arrays of scores and where they came from (i.e. "above", "left", or "diagonal")
            indMax = argmax(scoreArr)
            indMax==1 && push!(scoreDict, (i, j)=>[maxScore, "above"])
            indMax==2 && push!(scoreDict, (i, j)=>[maxScore, "left"])
            indMax==3 && push!(scoreDict, (i, j)=>[maxScore, "diagonal"])
        end
    end
    return scoreDict
end

"""
    function nwalign(s1, s2; match=1, mismatch=-1, gap=-1)

Takes two biological sequences as positional arguments, and a match, mismatch, and gap score as kwargs.
Uses the Needleman-Wunsch algorithm and returns one of the best global alignments of the two sequences as Tuple.
Gaps in the sequences are denoted with -'s.

Examples
≡≡≡≡≡≡≡≡≡
    julia> nwalign("AATTGGCC", "AAGGTTCC", mismatch=-2)
    ("AA--TTGGCC", "AAGGTT--CC")

    julia> nwalign("AATTGGCC", "AAGGTTCC", gap=-2)
    ("AATTGGCC", "AAGGTTCC")
"""
function nwalign(s1, s2; match=1, mismatch=-1, gap=-1)
    scoredMatDict = nwscorematrix(s1, s2; match=match, mismatch=mismatch, gap=gap)

    alignedSeq1 = ""
    alignedSeq2 = ""

    #traverse through matrix while i>1 and j>1
    #but first start at the cell in the bottom right corner
    i = length(s1) + 1
    j = length(s2) + 1
    while i>1 && j>1
        key = i, j
        direc = scoredMatDict[i, j][2]
        
        direc == "above" && (alignedSeq1 *= s1[i-1]; alignedSeq2 *= "-"; i-=1)
        direc == "left" && (alignedSeq1 *= "-"; alignedSeq2 *= s2[j-1]; j-=1)
        direc == "diagonal" && (alignedSeq1 *= s1[i-1]; alignedSeq2 *= s2[j-1]; i-=1; j-=1)
    end
    return reverse(alignedSeq1), reverse(alignedSeq2)
end

"""
    function alignmentHeatmap(s1, s2)
    
Takes two biological sequences as Strings and returns a matrix of size length(s1) by length(s2)
where cells corresponding to identical residues contain the value 1.0 and all other cells contain the value 0.0
Assumes that s1 is "lined" vertically, with s1[1] at matrix pos (1, 1) and s1[length(s1)] at matrix pos (length(s1), 1).
Assumes that s2 is "lined" horizontally, with s2[1] at matrix pos (1, 1) and s2[length(s2)] at matrix pos (1, length(s2)).

Useful for making a Bioinformatics dot-plot with Plots.jl (i.e. a monochromatic B&W heatmap
where 0.0 corresponds to white and 1.0 to black), as this function generates the data for
coloring each grid.

Examples
≡≡≡≡≡≡≡≡≡
    julia> alignmentHeatmap("AA", "GGA")[1, 3]
    1.0

    julia> alignmentHeatmap("AGTTTTAA", "AAGTGTGAC")
    8×9 Matrix{Float64}:
    1.0  1.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0
    0.0  0.0  1.0  0.0  1.0  0.0  1.0  0.0  0.0
    0.0  0.0  0.0  1.0  0.0  1.0  0.0  0.0  0.0
    0.0  0.0  0.0  1.0  0.0  1.0  0.0  0.0  0.0
    0.0  0.0  0.0  1.0  0.0  1.0  0.0  0.0  0.0
    0.0  0.0  0.0  1.0  0.0  1.0  0.0  0.0  0.0
    1.0  1.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0
    1.0  1.0  0.0  0.0  0.0  0.0  0.0  1.0  0.0
"""
function alignmentHeatmap(s1, s2)
    matrix = zeros(length(s1), length(s2))
    # if it's the same, color it => make the number in the cell a 1
    for i in 1:length(s1)
        for j in 1:length(s2)
            s1[i]==s2[j] && (matrix[i, j] = 1.0)
        end
    end
    return matrix
end
"""
    function mis_matchSeq(alignedSeq)

Takes a tuple of biological sequences already aligned. Assumes that the aligned sequences have the same length.
Returns a dictionary of 3 key-value pairs, where the keys are the Strings “Match#”, “Mismatch#”, “Mis/Match%” and 
corresponding values are an Int of the number of total positions with matching residues, 
an Int of the number of total positions with mismatching residues, and a Tuple containing the % of matches and the % of mismatches, in this order. 

Examples
≡≡≡≡≡≡≡≡≡
    julia> mis_matchSeq(("AGCT", "AGGT"))
    Dict{String, Any} with 3 entries:
    "Mis/Match%" => (0.75, 0.25)
    "Mismatch#"  => 1
    "Match#"     => 3

    julia> mis_matchSeq(("MFVFLVLLPLVS", "MFIFL-LFLTLT"))["Mismatch#"]
    7
"""
function mis_matchSeq(alignedSeq)
    s1 = alignedSeq[1]
    s2 = alignedSeq[2]
    matchNum = 0
    mismatchNum = 0

    for index in 1:length(s1) #assumes the 2 sequences have the same length
        if s1[index] == s2[index]
            matchNum += 1
        else
            mismatchNum += 1
        end
    end

    return Dict("Match#"=> matchNum, "Mismatch#"=> mismatchNum, "Mis/Match%" => (matchNum/length(s1), mismatchNum/length(s1)))
end 
"""
    function mis_matchSubsequences(alignedSeq)

Takes a tuple of biological sequences already aligned. Assumes that the aligned sequences have the same length.
Used to sift through two aligned sequences and return data on matching and mismatching subsequences.
Matching subsequences: consecutive, identical residues at the same positions in sequence 1 and 2
Mismatching subsequences: consecutive, non-identical residues at the same positions in sequence 1 and 2

Returns an Array of items of self-defined type AlignmentSubsequences with fields `type`, `positions`, and `sequences`.
The field `type` specifies whether the subsequences are a `"match"` or `"mismatch"`,
the field `positions` specifies the positions of the residues in the subsequence as a UnitRange,
the field `sequences` contains an Array of the matching subsequence or mismatching subsequences. 

Examples
≡≡≡≡≡≡≡≡≡
    julia> aligned = ("MFIFL-LFLTLTSGSAA", "MFVFLVL-LPLVS-SAA")
    julia> mis_matchSubsequences(aligned)
    13-element Vector{Any}:
    AlignmentSubsequences("match", 1:2, ["MF"])
    AlignmentSubsequences("mismatch", 3:3, ["I", "V"])
    AlignmentSubsequences("match", 4:5, ["FL"])
    ⋮
    AlignmentSubsequences("mismatch", 14:14, ["G", "-"])
    AlignmentSubsequences("match", 15:17, ["SAA"])

    julia> mis_matchSubsequences(("AGTCT", "AGACT"))[3].type
    "match"
    julia> mis_matchSubsequences(("AGTCT", "AGACT"))[3].positions
    4:5
    julia> mis_matchSubsequences(("AGTCT", "AGACT"))[3].sequences
    ["CT"]
"""
struct AlignmentSubsequences
    type::String # "match" or "mismatch"
    positions::UnitRange # e.g. 1:17
    sequences::Vector # e.g. ["AGCT"]
end
function mis_matchSubsequences(alignedSeq)
    arr = []
    s1 = alignedSeq[1]
    s2 = alignedSeq[2]
    matchStr = ""
    mismatchStr1 = ""
    mismatchStr2 = ""

    for index in 1:length(s1)
        s1[index] == s2[index] && (matchStr *= s1[index])
        s1[index] != s2[index] && (mismatchStr1 *= s1[index]; mismatchStr2 *= s2[index])

        
        if index < length(s1)
            # if current subsequences are matching and next ones are not, 
                # 1. make a struct
                # 2. push this struct to arr
                # 3. re-set matchStr to be an empty String 
            if s1[index] == s2[index] && s1[index+1] != s2[index+1]
                subseq = AlignmentSubsequences("match", index-length(matchStr)+1:index, [matchStr])
                #@info subseq
                push!(arr, subseq)
                matchStr = ""
            end

            # if current subsequences are mismatching and next ones are not, 
                # 1. make a struct
                # 2. push this struct to arr
                # 3. re-set mismatchStr1 and 2 to be empty Strings
            if s1[index] != s2[index] && s1[index+1] == s2[index+1]
                subseq = AlignmentSubsequences("mismatch", index-length(mismatchStr1)+1:index, [mismatchStr1, mismatchStr2])
                #@info subseq
                push!(arr, subseq)            
                mismatchStr1 = ""
                mismatchStr2 = ""
            end
        end
        
        if index == length(s1) # dealing with the very last subsequence
            if s1[index] == s2[index]
                subseq = AlignmentSubsequences("match", index-length(matchStr)+1:index, [matchStr])
                #@info subseq
                push!(arr, subseq)
            end

            if s1[index] != s2[index]
                subseq = AlignmentSubsequences("mismatch", index-length(mismatchStr1)+1:index, [mismatchStr1, mismatchStr2])
                #@info subseq
                push!(arr, subseq)
            end
        end
    end

    return arr
end



# Don't forget to export your functions!
end # module BioinformaticsBISC195


