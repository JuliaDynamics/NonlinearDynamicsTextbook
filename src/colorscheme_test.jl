# with help from 
# https://davidmathlogic.com/colorblind/

using PyPlot

COLORS = [
    # First color, it absolutely must be purple.
    "#6F4AC7", # most vibrant purple we can go for
    # "#745BB3", # between the above and below.
    # "#634e97", # original JuliaDynamics docs color
    # The next color is easy-peazy, because we cheat and use black!
    # The choice that remains is how light of a black hue we want
    # (pure black should be avoided if possible)
    # "#343434", # jet black
    "#000000", # pure black
    # "#2E2E2E",
    # One can only choose between Green or Orange to complement Purple.
    # I like Cyan more than Green, so I went for Orange as complement.
    # This establishes Cyan as second color:
    "#33CBD8",
    # Then, high brightness and saturation, red/orange hue. Possibilities:
    # "#DA3B3B",
    # "#E22411",
    "#E82727",
    # The last two colors are very tricky. One needs to be low saturation
    # and low brightness
    # Second last color needs to be dark and saturated. Possibilities:
    # "#8A4C75",
    # "#4C6C8A",
    # "#778CA0",
    "#6D7484",
    # "#84677A",
    # "#968A8A",
    # "#736C88",
    # "#806C88",
    # "#4C6C8A",
    # Last color is tricky. Needs to be bright, greenish. Possibilities:
    # "#2AC75D",
    # "#DCCE7D",
    # "#B6D840",
    "#A1D000",
]

if true # test color scheme
    fig = PyPlot.figure(figsize = (20, 15)) # show colors
    ax1 = subplot(231)
    ax2 = subplot(232)
    ax3 = subplot(233)
    ax4 = subplot(223)
    ax5 = subplot(224)
    lw = 60
    L = length(COLORS)
    for (i, c) in enumerate(COLORS)
        chsv = matplotlib.colors.rgb_to_hsv(matplotlib.colors.to_rgb(c))
        ax1.plot([0, 1], [0, 0] .+ i, color = c, lw = lw)
        ax1.set_title("color")
        ax2.plot([0, 1], [0, 0] .+ i, color = string(chsv[3]), lw = lw)
        ax2.set_title("brightness")
        ax3.plot([0, 1], [0, 0] .+ i, color = string(chsv[2]), lw = lw)
        ax3.set_title("saturation")
        x = 0:0.05:5π
        ax4.plot(x, cos.(x .+ i/2) .+ rand(length(x))/2; color=c, lw = 2)
        ax5.bar(collect(1:4) .+ (i-1)/L, 0.5rand(4) .+ 0.5, 1/L; color=c)
    end
    ax5.grid(false)
    ax4.grid(false)
    fig = tight_layout()
end
