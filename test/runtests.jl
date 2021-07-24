using BioinformaticsBISC195
using Test

@testset "BioinformaticsBISC195" begin
    
@testset "Using Strings" begin
    
    @testset "normalizeDNA" begin
        @test normalizeDNA("aatgn") == "AATGN"
        @test_throws Exception normalizeDNA("ZCA")
        @test_throws Exception normalizeDNA(42)
        c = normalizeDNA('C') 
        @test c == "C"
        @test typeof(c) == String
    end # normalizeDNA

    @testset "composition" begin
        seq = rand(['A','T','G','C','N'], 20) |> join
        bc = composition(seq)
        @test bc isa Tuple

        @test bc[1] == count(x-> x == 'A', seq)
        @test bc[3] == count(x-> x == 'C', seq)
        @test bc[2] == count(x-> x == 'G', seq)
        @test bc[4] == count(x-> x == 'T', seq)
        @test bc[5] == count(x-> x == 'N', seq)

        bc = composition(lowercase(seq))

        @test bc[1] == count(x-> x == 'A', seq)
        @test bc[3] == count(x-> x == 'C', seq)
        @test bc[2] == count(x-> x == 'G', seq)
        @test bc[4] == count(x-> x == 'T', seq)
        @test bc[5] == count(x-> x == 'N', seq)
    end # composition

    @testset "gc_content" begin
        @test gc_content("ANTG") == 0.25
        @test gc_content("cccggg") * 100 == 100.0
        @test gc_content("ATta") == 0.0
        @test_throws Exception gc_content("ATtz")
    end # gc_content

    @testset "complement" begin
        @test complement("ATTAN") == "TAATN"
        @test complement("gcta") == "CGAT"
        @test complement("nnnnnnn") == "NNNNNNN"
        @test_throws Exception complement("AZC")
    end # complement

    @testset "reverse_complement" begin
        @test reverse_complement("ATTAN") == "NTAAT"
        @test reverse_complement("gcta") == "TAGC"
        @test reverse_complement("nnnnnnn") == "NNNNNNN"
        @test_throws Exception reverse_complement("AZC")
    end # reverse_complement

    @testset "parse_fasta" begin
        testpath = normpath(joinpath(@__DIR__, "..", "data"))
        genomes = joinpath(testpath, "cov2_genomes.fasta")
        ex1_path = joinpath(testpath, "ex1.fasta")
        ex2_path = joinpath(testpath, "ex2.fasta")

        ex1 = parse_fasta(ex1_path)
        @test ex1 isa Tuple
        @test all(x-> x isa String, ex1[1])
        @test all(x-> x isa String, ex1[2])

        @test ex1[1] == ["ex1.1 | easy", "ex1.2 | multiline"]
        @test ex1[2] == ["AATTATAGC", "CGCCCCCCAGTCGGATT"]

        @test_throws Exception parse_fasta(ex2_path)

        cov2 = parse_fasta(genomes)
        @test length(cov2[1]) == 8
        @test length(cov2[2]) == 8
    end #parse_fasta

    @testset "removeUnder25k" begin
        parse1 = parse_fasta("../data/test1.fasta")
        parse2 = parse_fasta("../data/test2.fasta")
        cleaned1 = removeUnder25k(parse1)
        cleaned2 = removeUnder25k(parse2)
        @test length(findall(seq -> length(seq)<25000, cleaned2[2])) == 0
        @test length(cleaned1) ==length(parse1[2])
        @test length(cleaned2[1]) == length(cleaned2[2])
    end

    @testset "uniqueKmers" begin
        testex1 = uniqueKmers("AGTCTCTCTCGGGGTA", 4)
        @test length(testex1) == 10
        @test testex1 isa Vector
        @test testex1[end] isa String
        @test testex1[end] == "GGTA"
        @test length(testex1[3]) == 4

        testex2 = uniqueKmers("ATTTATATAGCTNNA", 3)
        @test length(testex2) == 8
        @test length(testex2[1]) == 3
        @test testex2[begin] == "ATT"
        @test testex2[3] == "TTA"
        @test testex2 isa Vector
    end

    @testset "compareKmerSets" begin
        testex1 = uniqueKmers("AGTCTCTCTCGGGGTA", 4)
        testex2 = uniqueKmers("ATTTATATAGCTNNA", 3)
        @test compareKmerSets(testex1, testex2) == 1.0
        @test compareKmerSets(testex1, testex2) isa Float64
        
        testex3 = ["ATT", "AGT", "CTC", "AGC"]
        testex4 = ["ATT", "CTC", "AGG", "AGA"]
        @test compareKmerSets(testex3, testex4) == 0.6666666666666667
    end

    @testset "nwscore" begin
        @test_throws ArgumentError nwscore(nothing, nothing)
        @test nwscore('A', 'G') == -1
        @test nwscore('A', 'A'; match = 10) == 10
        @test nwscore('A', 'T'; mismatch = -4) !=-1
        @test_throws MethodError nwscore("A", "C")
    end

    @testset "nwsetupmatrix" begin
        test5 = nwsetupmatrix("AGAGAGACT", "AGACAGCCT")
        @test test5 isa Matrix
        nrow, ncol = size(test5)
        @test nrow == ncol
        @test nrow == 10
        @test test5[1:10] isa Vector
        @test test5[1:10][10] == -9.0
        test6 = nwsetupmatrix("AGCT", "AGTTGAGAGTAT"; gap = -3)
        nrow2, ncol2 = size(test6)
        @test test6[1, 13] == -3 * (length("AGTTGAGAGTAT"))
        @test test6[5] == -12.0
        @test nrow2 == length("AGCT") + 1
        @test ncol2 == length("AGTTGAGAGTAT") + 1
    end

    @testset "nwscorematrix" begin
        ex7 = nwscorematrix("AGCT", "AGTTGAGAGTAT")
        @test ex7 isa Dict
        @test ex7[(4, 6)]==[-1.0, "left"]
        @test ex7[(2, 8)]==[-5.0, "left"]
        @test ex7[(3, 3)][1] isa Float64
        @test ex7[(3, 3)][2] isa String
        @test ex7[(5, 12)][2] == "left"
        @test collect(keys(ex7))[1] isa Tuple
        @test collect(keys(ex7))[2] == (5, 5)
        @test length(ex7) == 48
    end

    @testset "nwalign" begin
        ex8 = nwalign("AGCT", "AGCCA"; match = 2, mismatch = -3)
        @test ex8 isa Tuple
        @test ex8 == ("AGC--T", "AGCCA-")
        @test nwalign("AGCT", "AGCCA")[1] == "AGCT-"
        @test nwalign("AGCT", "AGCCA")[2] == "AGCCA" 

        ex9 = nwalign("GTGCCGA", "GAACCGTA"; gap=-2)
        @test ex9 == ("GTGCCG-A", "GAACCGTA")
    end

    @testset "alignmentHeatmap" begin
        ex10 = alignmentHeatmap("GTA", "ATA")
        nrow, ncol = size(ex10)
        @test nrow == ncol == 3
        @test ex10[3, 1] == 1.0
        @test ex10[2, 2] == 1.0
        @test ex10[3, 3] == 1.0

        ex11 = alignmentHeatmap("ATAAA", "ATAAA")
        @test ex11[1, 1] == ex11[2, 2] == ex11[3, 3] == ex11[4, 4] == ex11[5, 5]== 1.0
    end

    @testset "mis_matchSeq" begin
        aligned1 = nwalign("ATATAGCAAAAAA", "ATAGGCTATAAAA")
        ex12 = mis_matchSeq(aligned1)
        @test ex12 isa Dict
        @test ex12["Match#"] isa Int64
        @test ex12["Mismatch#"] isa Int64
        @test ex12["Mis/Match%"] isa Tuple
        @test ex12["Match#"] == 10
        @test ex12["Mismatch#"] == 4
        @test round(ex12["Mis/Match%"][1], digits=6) == 0.714286
    end

    @testset "mis_matchSubsequences" begin
        aligned1 = nwalign("ATATAGCAAAAAA", "ATAGGCTATAAAA")
        ex14 = mis_matchSubsequences(aligned1)
        @test ex14 isa Vector
        @test ex14[1].type == "match"
        @test length(ex14[1].sequences) == 1
        @test length(ex14[7].sequences[1]) == 4
        @test ex14[4].positions isa UnitRange
        @test ex14[4].positions == 8:8
        @test length(ex14[4].sequences) == 2
    end
end # strings

# @testset "Using BioSequences" begin
    
#     @testset "normalizeDNA" begin
#         @test normalizeDNA("aatgn") == dna"AATGN"
#         @test_throws Exception normalizeDNA("ZCA")
#         @test_throws Exception normalizeDNA(42)
#         c = normalizeDNA('C') 
#         @test c == dna"c"
#         @test c isa LongSequence
#     end #  normalizeDNA

#     @testset "gc_content" begin
#         @test gc_content(dna"ANTG") == 0.25
#         @test gc_content(dna"cccggg") * 100 == 100.0
#         @test gc_content(dna"ATta") == 0.0
#     end #  composition

#     @testset "composition" begin
#         seq = rand(['A','T','G','C','N'], 20) |> join |> LongDNASeq
#         bc = composition(seq)

#         @test bc[DNA_A] == count(==(DNA_A), collect(seq))
#         @test bc[DNA_C] == count(==(DNA_C), collect(seq))
#         @test bc[DNA_G] == count(==(DNA_G), collect(seq))
#         @test bc[DNA_T] == count(==(DNA_T), collect(seq))
#         @test bc[DNA_N] == count(==(DNA_N), collect(seq))
#     end #  gc_content

#     @testset "complement" begin
#         @test complement(dna"ATTAN") == dna"TAATN"
#         @test complement(dna"gcta") == dna"CGAT"
#         @test complement(dna"nnnnnnn") == dna"NNNNNNN"
#     end #  complement

#     @testset "reverse_complement" begin
#         @test reverse_complement(dna"ATTAN") == dna"NTAAT"
#         @test reverse_complement(dna"gcta") == dna"TAGC"
#         @test reverse_complement(dna"nnnnnnn") == dna"NNNNNNN"
#     end #  reverse_complement

#     @testset "parse_fasta" begin
#         testpath = normpath(joinpath(@__DIR__, "..", "data"))
#         genomes = joinpath(testpath, "cov2_genomes.fasta")
#         ex1_path = joinpath(testpath, "ex1.fasta")
#         ex2_path = joinpath(testpath, "ex2.fasta")

#         ex1 = parse_fasta(ex1_path)
#         @test ex1 isa Tuple
#         @test all(x-> x isa String, ex1[1])
#         @test all(x-> x isa LongSequence, ex1[2])

#         @test ex1[1] == ["ex1.1 | easy", "ex1.2 | multiline"]
#         @test ex1[2] == [dna"AATTATAGC", dna"CGCCCCCCAGTCGGATT"]

#         @test_throws Exception parse_fasta(ex2_path)

#         cov2 = parse_fasta(genomes)
#         @test length(cov2[1]) == 8
#         @test length(cov2[2]) == 8
#     end # parse_fasta

# end # BioSequences

end # BioinformaticsBISC195
