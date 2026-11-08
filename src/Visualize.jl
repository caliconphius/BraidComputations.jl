_svg_escape(x) = replace(string(x), "&"=>"&amp;", "<"=>"&lt;", ">"=>"&gt;", "\""=>"&quot;")

function _svg_document(body::AbstractString, width::Integer, height::Integer;
                        title="", orientation=:left_to_right)
    orientation in (:left_to_right, :top_to_bottom) ||
        throw(ArgumentError("orientation must be :left_to_right or :top_to_bottom"))
    if orientation == :top_to_bottom
        body = "<g transform=\"translate($height,0) rotate(90)\">$body</g>"
        width, height = height, width
    end
    title_markup = isempty(title) ? "" : "<title>$(_svg_escape(title))</title>"
    """<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="$width" height="$height" viewBox="0 0 $width $height" role="img">
$title_markup
<body style="font-family:sans-serif;background:white;fill:black">$body</body>
</svg>
"""
end

_svg_line(x1, y1, x2, y2; width=2, color="#202020") =
    "<line x1=\"$x1\" y1=\"$y1\" x2=\"$x2\" y2=\"$y2\" stroke=\"$color\" stroke-width=\"$width\"/>"

_svg_text(x, y, text; size=12, anchor="middle") =
    "<text x=\"$x\" y=\"$y\" font-size=\"$size\" text-anchor=\"$anchor\">$(_svg_escape(text))</text>"

_svg_path(d; color="#202020", width=2) =
    "<path d=\"$d\" fill=\"none\" stroke=\"$color\" stroke-width=\"$width\" stroke-linecap=\"round\"/>"

function _crossing_paths(x0, x1, y0, y1, sign)
    xm = (x0 + x1) / 2
    a = "M $x0 $y0 C $xm $y0 $xm $y1 $x1 $y1"
    b = "M $x0 $y1 C $xm $y1 $xm $y0 $x1 $y0"
    body = _svg_path(a; width=2) * _svg_path(b; width=6, color="white") * _svg_path(b)
    sign > 0 ? body : _svg_path(b; width=2) * _svg_path(a; width=6, color="white") * _svg_path(a)
end

"""
    braid_svg(x; orientation=:left_to_right)

Render a braid as a scalar-vector image. Strands are drawn left-to-right by
default and top-to-bottom when requested. Positive and inverse crossings are
indicated by continuous/broken overpasses and labelled Artin generators.
"""
function braid_svg(x::Braid; orientation=:left_to_right,
                   crossing_width=58, strand_spacing=32)
    crossing_width > 20 && strand_spacing >= 16 ||
        throw(ArgumentError("drawing dimensions are too small"))
    word = letters(x)
    n = nstrands(x)
    height = 24 + n * strand_spacing + 40
    width = 40 + max(length(word), 1) * crossing_width + 40
    left = 20
    body = ""

    for (column, generator) in enumerate(word == Int[] ? [0] : word)
        x0 = left + (column-1)*crossing_width
        x1 = x0 + crossing_width
        if generator == 0
            for strand in 1:n
                y = 12 + strand*strand_spacing
                body *= _svg_line(x0, y, x1, y)
            end
            continue
        end
        i = abs(generator)
        for strand in 1:n
            y = 12 + strand*strand_spacing
            if strand == i
                continue
            elseif strand == i+1
                continue
            end
            body *= _svg_line(x0, y, x1, y)
        end
        y1 = 12 + i*strand_spacing
        y2 = y1 + strand_spacing
        body *= _crossing_paths(x0, x1, y1, y2, generator)
        label = generator > 0 ? "σ$(i)" : "σ$(i)⁻"
        body *= _svg_text((x0+x1)/2, height-10, label)
    end
    for strand in 1:n
        body *= _svg_text(8, 17 + strand*strand_spacing, strand; anchor="end", size=11)
    end
    _svg_document(body, width, height; title="Braid on $n strands", orientation)
end

braid_svg(word::AbstractVector{<:Integer}, n::Integer; kwargs...) =
    braid_svg(BraidGroup(n)(word); kwargs...)

"""Render standard parabolic curves as nested ellipses around ordered punctures."""
function crs_svg(system::ParabolicSystem; n=nstrands(system),
                 orientation=:left_to_right, puncture_spacing=44,
                 singleton_punctures=Int[])
    n >= 2 || throw(ArgumentError("specify n ≥ 2 to draw punctures"))
    all(P -> nstrands(P) == n, system.curves) ||
        throw(ArgumentError("curve and disk strand counts do not agree"))
    width = 80 + n*puncture_spacing
    height = 180
    center_y = 86
    body = "<rect x=\"10\" y=\"28\" width=\"$(width-20)\" height=\"116\" rx=\"55\" fill=\"white\" stroke=\"#202020\" stroke-width=\"2\"/>"
    xs = collect(range(40, width-40; length=n))
    body *= _svg_line(xs[1], center_y, xs[end], center_y; width=1, color="#777777")

    for P in system.curves
        firstpoint, lastpoint = extrema(punctures(P))
        center = (xs[firstpoint] + xs[lastpoint])/2
        rx = (lastpoint-firstpoint)*puncture_spacing/2 + 16
        ry = 16 + 4*puncture_count(P)
        body *= "<ellipse cx=\"$center\" cy=\"$center_y\" rx=\"$rx\" ry=\"$ry\" fill=\"none\" stroke=\"#2a6f97\" stroke-width=\"2\"/>"
        body *= _svg_text(center, center_y-max(22, ry+5), "A[$(P.first),$(P.last)]"; size=11)
    end
    for i in singleton_punctures
        1 <= i <= n || throw(ArgumentError("singleton punctures must belong to 1:n"))
        xpos = xs[i]
        body *= "<ellipse cx=\"$xpos\" cy=\"$center_y\" rx=\"10\" ry=\"18\" fill=\"none\" stroke=\"#666666\" stroke-width=\"1.5\"/>"
    end
    for (i, xpos) in enumerate(xs)
        body *= "<circle cx=\"$xpos\" cy=\"$center_y\" r=\"4\" fill=\"#202020\"/>"
        body *= _svg_text(xpos, center_y+25, i; size=11)
    end
    isempty(system) && (body *= _svg_text(width/2, 164, "empty CRS (all punctures are outer loops)"; size=11))
    _svg_document(body, width, height; title="Standard canonical reduction system in B$n", orientation)
end

crs_svg(x::Braid; kwargs...) = crs_svg(canonical_reduction_system(x).standard; n=nstrands(x), kwargs...)

"""Draw the standard outer CRS, including one loop around each uncovered puncture."""
function outer_crs_svg(x::Braid; kwargs...)
    system = canonical_reduction_system(x).standard
    outers = filter(C -> !any(P -> C != P && issubset(C, P), system.curves), system.curves)
    covered = reduce(vcat, punctures.(outers); init=Int[])
    uncovered = setdiff(1:nstrands(x), covered)
    system = ParabolicSystem(outers)
    crs_svg(system; n=nstrands(x), singleton_punctures=uncovered, kwargs...)
end

function save_braid_svg(path::AbstractString, x::Braid; kwargs...)
    open(path, "w") do io
        write(io, braid_svg(x; kwargs...))
    end
    path
end

function save_crs_svg(path::AbstractString, x::Union{Braid,ParabolicSystem}; kwargs...)
    open(path, "w") do io
        write(io, crs_svg(x; kwargs...))
    end
    path
end

function save_outer_crs_svg(path::AbstractString, x::Braid; kwargs...)
    open(path, "w") do io
        write(io, outer_crs_svg(x; kwargs...))
    end
    path
end
