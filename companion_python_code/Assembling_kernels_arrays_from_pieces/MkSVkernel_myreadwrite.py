########################
# 06/29/2020 A. Alekseenko
# This module contains functions and
# subroutines that read and write data to hard drive
########################


############################
# my_read_svkernel (path)
#
# This subroutine reads an array
#  containing values of kernel based on a singular value discretization.
# this array is computed by a fortran code. This subroutine reads the array
#
##############################
def my_read_svkernel (path):
    #open binary file
    savefile=open(path, 'rb')

    #the first write is to record the first and laste index of the components of the SVKrnl array
    # the time variable is double precision and each write in fortran ads 32 bits before and after the written staff.
    # the 32 bits contain the number of bytes in the record.
    # so 4 bytes+4 bytes+ 4 bytes + 4bytes is 16
    bytes=savefile.read(16)
    # now that we got the first piece of information, let us unpack it using struct
    import struct
    data=struct.unpack('=llll',bytes)
    scr1=data[0]
    si_start=data[1]
    si_end = data[2]
    scr2=data[3]
    # the second piece of information and third long integers are the first single index and the last single index of the
    # recorded components. This gives us a way to calculuate the length of the singular value kernel array,
    # 4 bytes padding, 4 bytes the files size and 4 bytes
    # next we need to read fsize double reals into a sequence (not sure what it will be, list or tuple
    bytes1 = savefile.read(4)
    array_len=struct.unpack('=l',bytes1)  # this should be the length of the singular vector collision kernel array
    num_rec=array_len[0]/8
    # a quick consistency check
    if (si_end-si_start+1) != (num_rec):
        print("the numbers are inconsistent. Stop")
        exit()
    bytes = savefile.read(int(num_rec)*8)
    fmt='='+str(si_end-si_start+1)+"d"
    solscr=struct.unpack(fmt,bytes)
    # next we convert the sequence into a numpy array to return it.
    import numpy
    sol=numpy.array(solscr).reshape(si_end-si_start+1)
    # debug
    #sh=sol.shape
    # end debug
    # close the file
    savefile.close()
    return sol, (si_end-si_start+1)


############################
# my_write_svkernel (sol,path)
#
# This subroutine writes the  array sol on the hard drive. It
#  contains values of kernel based on a singular value discretization.
# this array is computed by a fortran code. The subroutine my_read_svkernel is used to read peices.
# of the array. THen the pieces are combined into a single array. Then this subroutine is called to write it on the hard drive
# in the same format
#
##############################
def my_write_svkernel (sol, path):
    #open binary file
    savefile=open(path, 'wb')
    # the first write is to record the first and last indices of the components of the SVKrnl array
    # the time variable is double precision and each write in fortran ads 32 bits before and after the written staff.
    # the 32 bits contain the number of bytes in the record.
    # so 4 bytes+4 bytes+ 4 bytes + 4bytes is 16
    import struct
    bytes = struct.pack('=4l', 8, 1, sol.shape[0], 8)  # packing the two long integer -- array sizes
    savefile.write(bytes)
    # next we need to record the SVKrnl. First we will need to figure the length of the record.  read fsize double reals into a sequence (not sure what it will be, list or tuple
    reclen = sol.shape[0] * 8  # the number of bytes used by the entire array
    bytes = struct.pack("=l", reclen)
    savefile.write(bytes)
    # now we write the actual array SVKrnl
    #### Next write the actual array records
    import numpy as np
    bytes = sol.tostring(order='F')
    savefile.write(bytes)
    ###############################
    # conclusing record size record
    reclen = sol.shape[0] * 8  # the number of bytes used by the entire array
    bytes = struct.pack("=l", reclen)
    savefile.write(bytes)
    ################################
    savefile.close()
    return
