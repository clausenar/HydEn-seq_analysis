#!/usr/bin/env python3

from fpdf import FPDF
pdf = FPDF()
import glob

#imagelist=[i.strip() for i in open('files.txt')]

imagelist=[i for i in glob.glob('./img/*.png')]

#imagelist=imagelist[0:4]

print (imagelist)
print (len(imagelist))
# imagelist is the list with all image filenames


for i in imagelist:
    print (i)


i=1

y=0

for image in imagelist:
    #print (image)
    strain_start=image.find('ARC')
    strain=image[strain_start:strain_start+6]+str('\n')
    print (strain)

    counter=i%7

    if i%7==1:
        print ("adding new page")
        pdf.add_page()
        y=0
        print ("new identier")
        pdf.set_font('Arial', '', 6)
        pdf.write(12,strain)
        y=60

    if i%7==2:
        y=90
    if i%7==3:
        y=120
    if i%7==4:
        y=150
    if i%7==5:
        y=180
    if i%7==6:
        y=210
    if i%7==0:
        y=240

    #pdf.write(3,image)
    pdf.write(3,image)
    pdf.image(image,0,y,90,30)


    print(i)
    print ("i//7",i//7)
    print ("i%7",i%7)
    print("")
    i+=1





pdf.output("collected_figures.pdf", "F")

