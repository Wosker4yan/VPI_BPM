#TAPERS, WAVEGUIDES, SBENDS
from diffractio import degrees, plt, um,mm


def waveguide(x, z,width ,length , effective_index, rotation, scalar_mask):
    scalar_mask.rectangle(
        r0=(x * um, z * um),
        size=(width * um, 2 * length * um),
        angle=rotation * degrees,
        refraction_index=effective_index)


def grating_coupler(x, z,period ,width, teeth_number,length_base, ff, etch_depth,index_teeth, effective_index, substrate_width, substrate_index, scalar_mask):
    length_full = teeth_number*period
    scalar_mask.rectangle(
        r0=(x * um, z * um),
        size=(width * um, 2 * length_base * um),
        angle=0 * degrees,
        refraction_index=effective_index)
    W =period*ff



    for i in range(teeth_number):
        scalar_mask.rectangle(
            r0=(0.5*width  * um, length_base + i*2*period * um),
            size=(etch_depth * um, 2 * W * um),
            angle=0 * degrees,
            refraction_index=index_teeth)

    scalar_mask.rectangle(
        r0=(x * um, 0+length_base * um),
        size=(width * um, 10*2*length_full*um),
        angle=0 * degrees,
        refraction_index=effective_index)

    scalar_mask.rectangle(
        r0=(0-0.5*width * um-0.5*substrate_width*um, z*um),
        size=(substrate_width * um, 10*length_base + 4*length_full*um),
        angle=0 * degrees,
        refraction_index=substrate_index)

    scalar_mask.draw_refraction_index(draw_borders=False)


