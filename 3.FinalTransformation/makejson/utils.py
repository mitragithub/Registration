import numpy
import glob

datadir = 'toBrian'

def readaffine(brainno):
	dat = numpy.loadtxt('%s/%s/%s_globalaffinetrans.txt' % (datadir,brainno,brainno))
	Rt = dat.reshape(4,3)
	R = Rt[:3,:]
	T = Rt[3,:]
	return (R,T)


def readrigid(brainno,sliceno):
	rtfile = glob.glob('%s/%s/%s_AAV_registered_rigidtrans/*_%.4d_rigidtrans.txt' % (datadir,brainno,brainno,int(sliceno)))
	if len(rtfile)>0:
		dat = numpy.loadtxt(rtfile[0]) #.replace('&','\&'))
		Rt = dat.reshape(3,2)
		R = Rt[:2,:]
		T = Rt[2,:]
		return (R,T)
	return (numpy.eye(2),[0,0])

if __name__=="__main__":
	R,T=readrigid('PMD1080',261)
	print R
	print T

	R,T=readaffine('PMD1080')
	print R
	print T


