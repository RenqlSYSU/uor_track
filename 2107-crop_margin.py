'''
these are two attempts to remove white edges from an image
but all have failed
'''

from PIL import Image
from skimage import io

def corp_margin(img):
    img2=img.sum(axis=2)
    (row,col)=img2.shape
    row_top=0
    raw_down=0
    col_top=0
    col_down=0
    for r in range(0,row):
        if img2.sum(axis=1)[r]<700*col:
            row_top=r
            break
        
    for r in range(row-1,0,-1):
        if img2.sum(axis=1)[r]<700*col:
            raw_down=r
            break
            
    for c in range(0,col):
        if img2.sum(axis=0)[c]<700*row:
            col_top=c
            break

    for c in range(col-1,0,-1):
        if img2.sum(axis=0)[c]<700*row:
            col_down=c
            break
                                                                                    
    new_img=img[row_top:raw_down+1,col_top:col_down+1,0:3]
    return new_img
  
#subprocess.run('convert /home/users/qd201969/ERA5-1HR/stat'+str(nv)+'.png -bordercolor white -trim',shell=True) 
#im = io.imread("/home/users/qd201969/ERA5-1HR/stat"+str(nv)+".png")
#img_re = corp_margin(im)
#io.imsave("/home/users/qd201969/ERA5-1HR/stat"+str(nv)+".png",img_re)

#img = Image.open("/home/users/qd201969/ERA5-1HR/stat"+str(nv)+".png")
#img_c = img.crop(img.getbbox())
#img_c.save("/home/users/qd201969/ERA5-1HR/crop_stat"+str(nv)+".png")
