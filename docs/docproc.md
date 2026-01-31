# Document processing

## Strip SVG outputs from a Literate-generated Jupyter notebook

- Use `JSON.jl` to read the `ipynb` file and remove the output data.
- Use `Base.format_bytes()` to show human-readable file size.

```julia
function strip_svg(nbpath)
    oldfilesize = filesize(nbpath)
    nb = open(JSON.parse, nbpath, "r")
    for cell in nb["cells"]
        !haskey(cell, "outputs") && continue
        for output in cell["outputs"]
            !haskey(output, "data") && continue
            datadict = output["data"]
            if haskey(datadict, "image/png") || haskey(datadict, "image/jpeg")
                delete!(datadict, "text/html")
                delete!(datadict, "image/svg+xml")
            end
        end
    end
    rm(nbpath; force=true)
    write(nbpath, JSON.json(nb, 1))
    @info "Stripped SVG in $(nbpath). The original size is $(Base.format_bytes(oldfilesize)). The new size is $(Base.format_bytes(filesize(nbpath)))."
    return nbpath
end
```

## Convert a Jupyter notebook into a Literate notebook

Adapted from [NBInclude.jl](https://github.com/JuliaInterop/NBInclude.jl). Separators `#---` were added for code blocks.

```julia
function to_literate(nbpath; shell_or_help=r"^\s*[;?]")
    nb = open(JSON.parse, nbpath, "r")
    jlpath = splitext(nbpath)[1] * ".jl"
    open(jlpath, "w") do io
        separator = ""
        for cell in nb["cells"]
            if cell["cell_type"] == "code"
                s = join(cell["source"])
                isempty(strip(s)) && continue # Jupyter doesn't number empty cells
                occursin(shell_or_help, s) && continue  # Skip cells with shell and help commands
                print(io, separator, "#---\n", s)  # Literate code block mark
                separator = "\n\n"
            elseif cell["cell_type"] == "markdown"
                txt = join(cell["source"])
                print(io, separator, "#===\n", txt, "\n===#")
                separator = "\n\n"
            end
        end
    end
    return jlpath
end
```

## Force display format to PNG

In the VSCode plot panel and [Literate.jl](https://github.com/fredrikekre/Literate.jl) notebooks, PNG images are generally smaller than SVG ones. To force plots to be shown as PNG images, you can use [DisplayAs.jl](https://github.com/tkf/DisplayAs.jl) to show objects in a chosen MIME type.

```julia
import DisplayAs.PNG
using Plots
plot(rand(6)) |> PNG
```

If you don't want to install another package, you could just use `display()`.

```julia
using Plots
PNG(img) = display("image/png", img)
plot(rand(6)) |> PNG
```
