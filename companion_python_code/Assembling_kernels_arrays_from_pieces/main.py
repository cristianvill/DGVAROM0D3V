########################
# 06/29/2020 A. Alekseenko
# This is code that manipulates with pre-computed
#  components of the singular vecot collision kernel
########################


import numpy as np
import MkSVkernel_myreadwrite as MkSVKrw


svker0, size0 = MkSVKrw.my_read_svkernel("F:\\temp\\cleaned100_MkKrnl_AddMxwl\\M41Trim0\\ch000\\AddMwxlCln100M41Tr0SVKrnl_ch000.dat")
svker0 = np.reshape(svker0, (size0,1))
MkSVKrw.my_write_svkernel(svker0,"AddMwxlCln100M41Tr0SVKrnl_k17.dat")


svker1, size1 = MkSVKrw.my_read_svkernel("F:\\temp\\cleaned100_MkKrnl_AddMxwl\\M41Trim0\\ch001\\AddMwxlCln100M41Tr0SVKrnl_ch001.dat")
svker1 = np.reshape(svker1, (size1,1))

svkerAll =  np.concatenate((svker0,svker1),axis=0)
MkSVKrw.my_write_svkernel(svkerAll,"AddMwxlCln100M41Tr0SVKrnl_k22.dat")

svker2, size2 = MkSVKrw.my_read_svkernel("F:\\temp\\cleaned100_MkKrnl_AddMxwl\\M41Trim0\\ch002\\AddMwxlCln100M41Tr0SVKrnl_ch002.dat")
svker2 = np.reshape(svker2, (size2,1))

svkerAll =  np.concatenate((svkerAll,svker2),axis=0)
MkSVKrw.my_write_svkernel(svkerAll,"AddMwxlCln100M41Tr0SVKrnl_k25.dat")


svker3, size3 = MkSVKrw.my_read_svkernel("F:\\temp\\cleaned100_MkKrnl_AddMxwl\\M41Trim0\\ch003\\AddMwxlCln100M41Tr0SVKrnl_ch003.dat")
svker3 = np.reshape(svker3, (size3,1))

svkerAll =  np.concatenate((svkerAll,svker3),axis=0)
MkSVKrw.my_write_svkernel(svkerAll,"AddMwxlCln100M41Tr0SVKrnl_k28.dat")


svker4, size4 = MkSVKrw.my_read_svkernel("F:\\temp\\cleaned100_MkKrnl_AddMxwl\\M41Trim0\\ch004\\AddMwxlCln100M41Tr0SVKrnl_ch004.dat")
svker4 = np.reshape(svker4, (size4,1))

svkerAll =  np.concatenate((svkerAll,svker4),axis=0)
MkSVKrw.my_write_svkernel(svkerAll,"AddMwxlCln100M41Tr0SVKrnl_k30.dat")


svker5, size5 = MkSVKrw.my_read_svkernel("F:\\temp\\cleaned100_MkKrnl_AddMxwl\\M41Trim0\\ch005\\AddMwxlCln100M41Tr0SVKrnl_ch005.dat")
svker5 = np.reshape(svker5, (size5,1))

svkerAll =  np.concatenate((svkerAll,svker5),axis=0)

svker5, size5 = MkSVKrw.my_read_svkernel("F:\\temp\\cleaned100_MkKrnl_AddMxwl\\M41Trim0\\ch006\\AddMwxlCln100M41Tr0SVKrnl_ch006.dat")
svker5 = np.reshape(svker5, (size5,1))

svkerAll =  np.concatenate((svkerAll,svker5),axis=0)

svker5, size5 = MkSVKrw.my_read_svkernel("F:\\temp\\cleaned100_MkKrnl_AddMxwl\\M41Trim0\\ch007\\AddMwxlCln100M41Tr0SVKrnl_ch007.dat")
svker5 = np.reshape(svker5, (size5,1))

svkerAll =  np.concatenate((svkerAll,svker5),axis=0)
MkSVKrw.my_write_svkernel(svkerAll,"AddMwxlCln100M41Tr0SVKrnl_k36.dat")

svker6, size6 = MkSVKrw.my_read_svkernel("F:\\temp\\cleaned100_MkKrnl_AddMxwl\\M41Trim0\\ch008\\AddMwxlCln100M41Tr0SVKrnl_ch008.dat")
svker6 = np.reshape(svker6, (size6,1))

svkerAll =  np.concatenate((svkerAll,svker6),axis=0)
MkSVKrw.my_write_svkernel(svkerAll,"AddMwxlCln100M41Tr0SVKrnl_k37.dat")

svker6, size6 = MkSVKrw.my_read_svkernel("F:\\temp\\cleaned100_MkKrnl_AddMxwl\\M41Trim0\\ch009\\AddMwxlCln100M41Tr0SVKrnl_ch009.dat")
svker6 = np.reshape(svker6, (size6,1))

svkerAll =  np.concatenate((svkerAll,svker6),axis=0)

svker6, size6 = MkSVKrw.my_read_svkernel("F:\\temp\\cleaned100_MkKrnl_AddMxwl\\M41Trim0\\ch010\\AddMwxlCln100M41Tr0SVKrnl_ch010.dat")
svker6 = np.reshape(svker6, (size6,1))

svkerAll =  np.concatenate((svkerAll,svker6),axis=0)

svker6, size6 = MkSVKrw.my_read_svkernel("F:\\temp\\cleaned100_MkKrnl_AddMxwl\\M41Trim0\\ch011\\AddMwxlCln100M41Tr0SVKrnl_ch011.dat")
svker6 = np.reshape(svker6, (size6,1))

svkerAll =  np.concatenate((svkerAll,svker6),axis=0)

svker6, size6 = MkSVKrw.my_read_svkernel("F:\\temp\\cleaned100_MkKrnl_AddMxwl\\M41Trim0\\ch012\\AddMwxlCln100M41Tr0SVKrnl_ch012.dat")
svker6 = np.reshape(svker6, (size6,1))

svkerAll =  np.concatenate((svkerAll,svker6),axis=0)

svker7, size7 = MkSVKrw.my_read_svkernel("F:\\temp\\cleaned100_MkKrnl_AddMxwl\\M41Trim0\\ch013\\AddMwxlCln100M41Tr0SVKrnl_ch013.dat")
svker7 = np.reshape(svker7, (size7,1))

svkerAll =  np.concatenate((svkerAll,svker7),axis=0)

svker7, size7 = MkSVKrw.my_read_svkernel("F:\\temp\\cleaned100_MkKrnl_AddMxwl\\M41Trim0\\ch014\\AddMwxlCln100M41Tr0SVKrnl_ch014.dat")
svker7 = np.reshape(svker7, (size7,1))

svkerAll =  np.concatenate((svkerAll,svker7),axis=0)

svker7, size7 = MkSVKrw.my_read_svkernel("F:\\temp\\cleaned100_MkKrnl_AddMxwl\\M41Trim0\\ch015\\AddMwxlCln100M41Tr0SVKrnl_ch015.dat")
svker7 = np.reshape(svker7, (size7,1))

svkerAll =  np.concatenate((svkerAll,svker7),axis=0)
MkSVKrw.my_write_svkernel(svkerAll,"AddMwxlCln100M41Tr0SVKrnl_k45.dat")


svker7, size7 = MkSVKrw.my_read_svkernel("F:\\temp\\cleaned100_MkKrnl_AddMxwl\\M41Trim0\\ch016\\AddMwxlCln100M41Tr0SVKrnl_ch016.dat")
svker7 = np.reshape(svker7, (size7,1))

svkerAll =  np.concatenate((svkerAll,svker7),axis=0)

svker7, size7 = MkSVKrw.my_read_svkernel("F:\\temp\\cleaned100_MkKrnl_AddMxwl\\M41Trim0\\ch017\\AddMwxlCln100M41Tr0SVKrnl_ch017.dat")
svker7 = np.reshape(svker7, (size7,1))

svkerAll =  np.concatenate((svkerAll,svker7),axis=0)

svker7, size7 = MkSVKrw.my_read_svkernel("F:\\temp\\cleaned100_MkKrnl_AddMxwl\\M41Trim0\\ch018\\AddMwxlCln100M41Tr0SVKrnl_ch018.dat")
svker7 = np.reshape(svker7, (size7,1))

svkerAll =  np.concatenate((svkerAll,svker7),axis=0)
MkSVKrw.my_write_svkernel(svkerAll,"AddMwxlCln100M41Tr0SVKrnl_k48.dat")

svker7, size7 = MkSVKrw.my_read_svkernel("F:\\temp\\cleaned100_MkKrnl_AddMxwl\\M41Trim0\\ch019\\AddMwxlCln100M41Tr0SVKrnl_ch019.dat")
svker7 = np.reshape(svker7, (size7,1))

svkerAll =  np.concatenate((svkerAll,svker7),axis=0)

MkSVKrw.my_write_svkernel(svkerAll,"AddMwxlCln100M41Tr0SVKrnl_kxx.dat")


print(range(5))