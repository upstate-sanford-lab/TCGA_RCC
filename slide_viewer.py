import openslide
import os

filedir='/data/TCGA_RCC_sorted'
dirname='TCGA-CZ-5460'

#non-clear appearance (TCGA-CZ-5460)

def slide_viewer():
    for file in os.listdir(os.path.join(filedir,dirname)):
        if file.endswith('.svs'):
            filepath=os.path.join(filedir,dirname,file)
    oslide = openslide.OpenSlide(filepath)
    thumb = oslide.get_thumbnail(size=(1000, 1000))
    thumb.show()

slide_viewer()