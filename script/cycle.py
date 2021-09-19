import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np

fig=plt.figure(figsize=(5,5),dpi=300)   #创建画布

artists = []
for j in range(3):
    axe, = plt.subplot(1,1,1)    #创建子图
    axe.set_xticks(np.arange(0,j,1))
    axe.set_yticks(np.arange(0,j,1))
    axe.set_title('test '+str(j))   #设置子图标题
    artists.append(axe)

ani = animation.ArtistAnimation(fig=fig, artists=artists, repeat=True, interval=100)
plt.show()
ani.save('./2.mp4', fps=30)
#controls the spacing between subgraphs, and hspace is the vertical spacing
#plt.subplots_adjust(wspace=0.3,hspace=0.4)
#fig.savefig('./test.png',dpi=300) #保存图片
