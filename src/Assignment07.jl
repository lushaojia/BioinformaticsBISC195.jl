module Assignment07

export normalizeDNA,
       composition,
       gc_content,
       complement,
       reverse_complement


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

# Your code here.
# Don't forget to export your functions!

end # module Assignment07

