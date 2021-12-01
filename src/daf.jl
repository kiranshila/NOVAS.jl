import Base.show

export DAFArray, DAFFile, read_daf

const RECORD_LENGTH = 1024
const WORD_SIZE = 8
const SRH_SIZE = 24

struct FileRecord
    locidw::String
    nd::Int
    ni::Int
    locifn::String
    fward::Int
    bward::Int
    free::Int
    locfmt::String
end

# Don't read past the lcfmt, the ftpstring kills the reader as it has nonascii
function FileRecord(f::FortranFile)
    return FileRecord(read(f, FString{8}, Int32, Int32, FString{60}, Int32, Int32, Int32,
                           FString{8}; rec=1)...)
end

struct Summary
    doubles::Vector{Float64}
    ints::Vector{Int32}
end

function summaries(f::FortranFile, fr::FileRecord, nsum::Int, record::Int)
    types = Base.Iterators.flatten([((Float64, fr.nd), (Int32, fr.ni)) for _ in 1:nsum])
    s = read(f, Float64, Float64, Float64, types...; rec=record)[4:end]
    return [Summary(s[i], s[i + 1]) for i in 1:2:(2nsum)]
end

struct DAFArray
    doubles::Vector{Float64}
    ints::Vector{Int32}
    data::Vector{Float64}
end

struct DAFFile
    name::String
    id::String
    arrays::Vector{DAFArray}
end

function show(io::IO, file::DAFFile)
    println(io, "A DAF File")
    println(io, "Name             : ", file.name)
    println(io, "ID               : ", file.id)
    return println(io, "Number of Arrays : ", length(file.arrays))
end

which_record(addr) = ((addr - 1) >> 7) + 1

"""
    read_daf("filename.daf")

Reads an NAIF Double Precision Array Files (DAF) into a Julia `DAFFile`
"""
function read_daf(filename)
    # Figure out endianness
    f = FortranFile(filename; convert="little-endian", access="direct", recl=RECORD_LENGTH)
    # Read header
    fr = FileRecord(f)
    if fr.locfmt == "BIG-IEEE"
        @info "DAF File is Big Endian"
        f = FortranFile(filename; convert="big-endian", access="direct", recl=RECORD_LENGTH)
        fr = FileRecord(f)
    end
    # Collect all summaries and names
    summs = Summary[]
    sum_ptr = fr.fward
    while true
        # Read header
        next, prev, nsum = read(f, Float64, Float64, Float64; rec=sum_ptr)
        append!(summs, summaries(f, fr, Int(nsum), sum_ptr))
        # Move the ptr
        if iszero(next)
            break
        else
            sum_ptr = Int(next)
        end
    end
    # Collect the data
    data = []
    for summary in summs
        start_idx = summary.ints[end - 1]
        end_idx = summary.ints[end]
        record_range = which_record(start_idx):which_record(end_idx)
        vect = Float64[]
        for record in record_range
            rec = read(f, (Float64, 128); rec=record)
            a, b = 1, 128
            if record == record_range[1]
                a = ((start_idx - 1) % 128) + 1
            end
            if record == record_range[end]
                b = ((end_idx - 1) % 128) + 1
            end
            append!(vect, rec[a:b])
        end
        push!(data, vect)
    end
    # Construct our thing
    daf_arrays = [DAFArray(summs[i].doubles, summs[i].ints[1:(end - 2)], data[i])
                  for i in 1:length(summs)]
    return DAFFile(fr.locifn, fr.locidw, daf_arrays)
end
