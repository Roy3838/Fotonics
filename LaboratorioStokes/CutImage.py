import cv2
import numpy as np

Images = [0,0,0,0,0,0]
Centeres = [0,0,0,0,0,0]
Mark = [0,0,0,0,0,0]

#Images centers in order
centers = [[379,220],[376,270],[384,241],[356,203],[371,213],[374,215]]
#Size of cutted image
size = [200,150]

for j in range(1,2):
    for i in range(1,7):
        print("Imagen"+str(j)+"_"+str(i))
        img = cv2.imread("Images/"+str(j)+"_"+str(i)+".jpg")
        
        rows,cols,_=img.shape
        print("Rows",rows)
        print("Cols",cols)

        center = centers[i-1]
        cutted = img[center[1]-size[1]:center[1]+size[1],center[0]-size[0]:center[0]+size[0]]

        markhv = cutted.copy()
        markhv = cv2.line(markhv, (size[0],0), (size[0],size[1]*2), (0,0,255), 1)
        markhv = cv2.line(markhv, (0,size[1]), (size[0]*2,size[1]), (0,0,255), 1)
        
        rows,cols,_=markhv.shape
        print("Rows",rows)
        print("Cols",cols)

        Images[i-1] = img
        Centeres[i-1] = cutted
        Mark[i-1] = markhv

        cv2.imshow("Marked"+str(j)+"_"+str(i),markhv)
        cv2.imshow("Cuted"+str(j)+"_"+str(i),cutted)
        cv2.imshow("image"+str(j)+"_"+str(i),img)

        cv2.imwrite("ImagesCut/"+str(j)+"_"+str(i)+".jpg", cutted)
        cv2.imwrite("ImagesCut/Centers/"+str(j)+"_"+str(i)+".jpg", markhv)
        cv2.waitKey(0)
        cv2.destroyAllWindows()




