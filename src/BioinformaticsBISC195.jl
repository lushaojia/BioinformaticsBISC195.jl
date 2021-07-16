module BioinformaticsBISC195

export normalizeDNA,
       composition,
       gc_content,
       complement,
       reverse_complement,
       parse_fasta

# # uncomment the following line if you intend to use BioSequences types
# using BioSequences

"""
    normalizeDNA(::AbstractString)

Ensures that a sequence only contains valid bases
(or `'N'` for unknown bases).
Returns a String.
"""
function normalizeDNA(seq)
    seq = uppercase(string(seq))
    for base in seq
        # note: `N` indicates an unknown base
        occursin(base, "AGCTN") || error("Invalid base, $base")
    end
    return seq # change to `return LongDNASeq(seq)` if you want to try to use BioSequences types
end

"""
    composition(sequence)

Counts the number of each type of base
in a DNA sequence and returns a dictionary with the keys 'A', 'C',
'G', 'T', and 'N' whose corresponding values are their counts.

Throws an error when an invalid base (not 'A', 'C',
'G', 'T', or 'N') is encountered.

Examples  
≡≡≡≡≡≡≡≡≡≡
    julia> composition("ACCGGGTTTTN")
        Dict{Char,Int64} with 5 entries:
            'A' => 1
            'G' => 3
            'T' => 4
            'N' => 1
            'C' => 2
    

    julia> composition("AAX")
    ERROR: Invalid base, X
"""
function composition(sequence)
    normalizedSeq = normalizeDNA(sequence)
    A = G = C = T = N = 0
    for base in normalizedSeq
        if base == 'A'
            A += 1
        elseif base == 'G'
            G += 1
        elseif base == 'C'
            C += 1
        elseif base == 'T'
            T += 1
        else
            N += 1
        end
    end
    compDict = Dict('A'=>A, 'G'=>G, 'T'=>T, 'C'=>C, 'N'=>N)

    return compDict

end

"""
    gc_content(seq)

Calculates the GC ratio (the total number of G and C bases 
divided by the total length of the sequence) of a DNA sequence.

Accomodates sequences with ambiguous bases (e.g. base N)

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
    gc_ratio = (baseCount['G']+baseCount['C'])/length(seq)
    return gc_ratio
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

Accomodates sequences with ambiguous bases (e.g. base N)

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

Accomodates sequences with ambiguous bases (e.g. base N)

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
    function parse_fasta(path)

Reads a fasta-formated file and returns 2 vectors,
one containing parsed headers 
(which do not contain the character '>' at the beginning),
the other containing the entire sequence as a `String`.

The function validates sequences to see if they contain DNA bases only. Accomodates 
sequences with ambiguous bases (e.g. bases R, Y, S, W, K, M, B, D, H, V, N and gaps . or -).

Example
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

function parse_fasta(path)
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
            for base in line
                occursin(base, "AGCTRYSWKMBDHVN.-") || error("Invalid base, $base")
            end
            seqStr *= line
        end
    end

    # push the last conactenated DNA sequence to seqArr
    push!(seqArr, seqStr)

    return headerArr, seqArr
end



# Don't forget to export your functions!

end # module BioinformaticsBISC195

