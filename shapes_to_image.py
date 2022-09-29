import matplotlib.pyplot as plt

from vpipda import units
from vpipda.material import *
from vpipda.layout import *
from vpipda.device import *
from PIL import Image


#WILL NOT ATTEMTP THIS
#SBEND figure out how to make a straigth with normal cropping
# Define the substrate material:
ns = Dielectric(n=1, color='white')
# Let's define the same core material with two
# different colors, to better visualize
# 3D device layouts:
nc2 = Dielectric(n=1, color='k')
# Define cross-section generator:
#UNJUSTIFIELBLY COMPLICATED ANYWAY WE JUST NEED AN IMAGE
def xs(width):
    # Waveguide layout:
    substrate = HalfPlane(ns)
    rib = Rectangle(nc2, w=width, h='0.5um')
    return Layout2D(substrate, rib)




dev = SBendCosine(xs, width=1, length=10, height=5, port_names=[])
dev.plot()
plt.axis('off')

plt.show()



wg = Straight(xs, width=1, length=10, port_names=[])
plt.tight_layout()
wg.plot()
plt.axis('off')

plt.show()
filename = 'mask'
plt.savefig(filename, bbox_inches='tight', pad_inches=0)

im = Image.open(filename)



im_resized = im.resize(5000, Image.ANTIALIAS)
im_resized.save(filename, 'png')

plt.close()


