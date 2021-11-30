using FortranFiles, DataDeps
import Base.show

register(DataDep("de440", "",
                 "https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de440.bsp",
                 "a4ce9bf9b3282becc9f4b2ac3cebe03a2ae7599981aabd7265fd8482fff7c4b5"))

const RECORD_LENGTH = 1024
const WORD_SIZE = 8
const SRH_SIZE = 24

struct FileRecord
    LOCIDW::NTuple{8,UInt8}
    ND::UInt32
    NI::UInt32
    LOCIFN::NTuple{60,UInt8}
    FWARD::UInt32
    BWARD::UInt32
    FREE::UInt32
    LOCFMT::NTuple{8,UInt8}
    PRENUL::NTuple{603,UInt8}
    FTPSTR::NTuple{28,UInt8}
    PSTNUL::NTuple{297,UInt8}
end

struct SummaryRecordHeader
    next::Float64
    prev::Float64
    count::Float64
end

struct Summary{D,I}
    doubles::NTuple{D,Float64}
    ints::NTuple{I,Int32}
end

struct DAFArray
    name::String
    doubles::Vector{Float64}
    ints::Vector{Int32}
    data::Vector{Float64}
end

struct DAFFile
    name::String
    id::String
    comments::String
    arrays::Vector{DAFArray}
end

function show(io::IO,file::DAFFile)
    println(io, "A DAF File")
    println(io, "Name             : ", file.name)
    println(io, "ID               : ", file.id)
    println(io, "Number of Arrays : ", file.arrays |> length)
end

read_record(f::FortranFile, r::Int) = read(f, (UInt8, RECORD_LENGTH); rec=r)

which_record(addr) = ((addr - 1) >> 7) + 1

function read_daf(filename)
    # DataDep is specifically little-endian
    f = FortranFile(filename; convert="little-endian", access="direct", recl=RECORD_LENGTH)
    fr = reinterpret(FileRecord, read_record(f, 1))[1]
    # Determine how many comment records
    cmt_idxs = 2:(fr.FWARD - 1)
    cmnts = ""
    # Collect comments
    for idx in cmt_idxs
        cmnts *= strip(replace(String([read_record(f, idx)...]), "\0" => "\n"))
    end
    # Collect all summaries and names
    summary_type = Summary{Int(fr.ND),2 * ((fr.NI + 1) ÷ 2)}
    summaries = summary_type[]
    nms = String[]
    sum_ptr = Int(fr.FWARD)
    while true
        # Read summary
        summary_record = read_record(f, sum_ptr)
        sum_header = reinterpret(SummaryRecordHeader, summary_record[1:SRH_SIZE])[1]
        summaries_type = NTuple{Int(sum_header.count),summary_type}
        append!(summaries,
                reinterpret(summaries_type,
                            summary_record[(SRH_SIZE + 1):(sizeof(summaries_type) + SRH_SIZE)])[1])
        # Read names
        name_record = read_record(f, sum_ptr + 1)
        names_type = NTuple{Int(sum_header.count),NTuple{sizeof(summary_type),UInt8}}
        names′ = reinterpret(names_type, name_record[1:sizeof(names_type)])[1]
        for name in names′
            push!(nms, strip(String(collect(name))))
        end
        # Move the ptr
        if iszero(sum_header.next)
            break
        else
            sum_ptr = sum_header.next
        end
    end
    # Collect the data
    data = []
    for summary in summaries
        start_idx = summary.ints[end - 1]
        end_idx = summary.ints[end]
        record_range = which_record(start_idx):which_record(end_idx)
        vect = Float64[]
        for record in record_range
            rec = reinterpret(Float64, read_record(f, record))
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
    daf_arrays = [DAFArray(nms[i], summaries[i].doubles |> collect, summaries[i].ints[1:(end - 2)] |> collect,
                           data[i]) for i in 1:length(summaries)]
    dafname = String(fr.LOCIFN |> collect) |> strip
    dafid = String(fr.LOCIDW |> collect) |> strip
    return DAFFile(dafname, dafid, cmnts, daf_arrays)
end

de440 = read_daf("/home/kiran/.julia/datadeps/de440/de440.bsp")