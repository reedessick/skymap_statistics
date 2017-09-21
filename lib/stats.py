usage="""a module to compute basic statistics of skymaps"""

import healpy as hp
import numpy as np

#=================================================
#
# general helper methods
#
#=================================================
def cos_dtheta(theta1, phi1, theta2, phi2, safe=False):
	"""
	computes the angular separation between two points
	support arrays assuming they all have the same shape
	if safe: 
		we check to make sure that numerical error doesn't put out outside the range of cosine
	"""
	cosdtheta = np.cos(theta1)*np.cos(theta2) + np.sin(theta1)*np.sin(theta2)*np.cos(phi1-phi2)

	if safe:
		if isinstance(cosdtheta, np.ndarray):
			cosdtheta[cosdtheta<-1] = -1
			cosdtheta[cosdtheta>1] = 1
		else:
			if cosdtheta < -1:
				cosdtheta = -1
			elif cosdtheta > 1:
				cosdtheta = 1

	return cosdtheta

#=================================================
#
# methods involving manipulating a posterior
#
#=================================================
def rankmap(posterior, npix=None, normed=True):
        """
        converts a posterior into a rank map. Small ranks correspond to large posterior weight
        WARNING: rankmaps do not sum to unity!
        """
	if not npix:
	        npix = len(posterior)
        rankmap = np.empty(npix,int)
        rankmap[posterior.argsort()] = np.arange(npix)[::-1]
        if normed:
            rankmap /= 1.0*npix
        return rankmap

###
def resample(posterior, new_nside, nest=False):
	"""
	creates a new posterior with len=hp.nside2npix(new_nside)
		if new_nside > nside: we assign each new pixel an equal fraction of the parent pixel's value
		if new_nside < nside: we assign each new pixel the sum of the contained pixels
	"""
	return hp.ud_grade(posterior, new_nside, power=-2)

###
def __to_cumulative(posterior):
	"""
	returns a map corresponding to cumulative probabilities at each pixel
	assumes ``greedy binning'' algorithm
	"""
        ranking = posterior.argsort()[::-1] ### find out ordering
        cum = np.empty_like(posterior)
        cum[ranking] = np.cumsum( posterior[ranking] ) ### and assign it to the correct indecies within cum
        return cum
 
        ### this looks really slow...
#	cum = np.zeros_like( posterior )
#	c = 0.0
#	for i in posterior.argsort()[::-1]:
#		c += posterior[i]
#		cum[i] = c
#	return cum

###
def credible_region(posterior, conf):
	"""
	returns a list of pixels that correspond the minimum credible region
	"""
	if isinstance(conf, (int,float)):
		conf = np.array([conf])
	elif not isinstance(conf, np.ndarray):
		conf = np.array(conf)
	if np.any(conf < 0) or np.any(conf > 1):
		raise ValueError("conf must be between 0 and 1")
	
	cum = __to_cumulative(posterior)
	
	ind = np.arange(len(posterior))
	return [ind[cum<=c] for c in conf] ### return indecies corresponding to confidence levels

###
def p_value(posterior, theta, phi, nside=None):
	"""
	computes the p-value at which a given point was found
	"""
	if not nside:
		nside = hp.npix2nside(len(posterior))
	return np.sum(posterior[posterior>=posterior[hp.ang2pix(nside, theta, phi)]])

###
def estang(posterior, nside=None):
	"""
	returns the position associated with the maximum of the posterior
	"""
	if not nside:
		nside = hp.npix2nside(len(posterior))
	return hp.pix2ang(nside, posterior.argmax())

###
def searched_area(posterior, theta, phi, nside=None, nest=False, degrees=False):
	"""
	computes the searched area given a location
	"""
	if not nside:
		nside = hp.npix2nside(len(posterior))
	ipix = hp.ang2pix(nside, theta, phi, nest=nest)
	return np.sum(posterior>=posterior[ipix])*hp.nside2pixarea(nside, degrees=degrees)

###
def est_cos_dtheta(posterior, theta, phi, safe=False):
	"""
	returns the angular separation between the maximum of the posterior and theta, phi
	"""
	t, p = estang(posterior)
	return cos_dtheta(theta, phi, t, p, safe=safe)

###
def min_cos_dtheta(posterior, theta, phi, nside=None, nest=False, safe=False):
	"""
	computes the maximum angular separation between any point in the area with p > p(theta, phi) and the estimated position
	"""
	if not nside:
		nside = hp.npix2nside(len(posterior))
	ipix = hp.ang2pix(nside, theta, phi, nest=nest)
	## get all ang included in this area
	thetas, phis = hp.pix2ang(nside, np.arange(len(posterior))[posterior >= posterior[ipix]])

	t, p = estang(posterior)

	return np.min(cos_dtheta(thetas, phis, t, p, safe=safe))

###
def min_all_cos_dtheta(pix, nside, nest=False, safe=False):
	"""
	computes the maximum angular separation between any two pixels within pix=[ipix,ipix,...]
	does this by direct comparison (N^2 computation) between pixel centers
	"""
	min_c = 1.0
	thetas, phis = hp.pix2ang(nside, pix)
	for i, (t1, p1) in enumerate(zip(thetas, phis)[:-1]):
		t2 = thetas[i+1:] 
		p2 = phis[i+1:]
		c = np.amin(cos_dtheta(t1,p1,t2,p2,safe=safe))
		if c < min_c:
			min_c = c
	return min_c

###
def min_all_cos_dtheta_fast(pix, nside, nest=False, safe=False):
	"""
	computes the maximum angular separation between any two pixels within pix=[ipix,ipix,...]
	does this with a boarder-to-boarder search after checking some other things
	this could be in error up to the pixel size (hopefully small)

	This algorithm is due in part to Antonios Kontos, who helped Reed Essick think through the details
	"""
	Npix = len(pix)
	### check to see if more than half the sky is filled
	if Npix*hp.nside2pixarea( nside ) > 2*np.pi:
		return -1 
	### check to see if the antipode of any point is in the set
	npix =  hp.nside2npix( nside )
	selected = np.zeros( npix, dtype=bool )
	selected[pix] = True
	antipodes = np.zeros_like( selected )
	vec = -np.array(hp.pix2vec( nside, pix, nest=nest )) ### reflect vectors
	antipodes[ hp.vec2pix( nside, vec[0], vec[1], vec[2], nest=nest ) ] = True ### reflection -> antipode
	                                                            ### reflection accomplished in cartesian coords
	if np.sum( selected*antipodes ): ### point and it's antipode are in the set
		return -1 ### could be in error by the pixel size...

	### boarder-to-boarder search
	boarder_pix = []
	for bpix in __into_boarders( nside, pix, nest=nest ):
		boarder_pix += bpix
	return min_all_cos_dtheta( boarder_pix, nside, nest=nest, safe=False )

def __into_boarders(nside, pix, nest=False):
        """
        extracts the boarder from the list of pixels (pix)
        """
        ### establish an array representing the included pixels
        npix = hp.nside2npix(nside)

        truth = np.zeros((npix,), bool)
        truth[pix] = True

        abstruth = np.zeros((npix,), bool) ### need completely separate object, not just a pointer
        abstruth[pix] = True

        pixnums = np.arange(npix) ### convenient array we establish once

        boarders = []
        while truth.any():
                ipix = pixnums[truth][0] ### take the first pixel
                truth[ipix] = False ### remove it from the global set
                boarder = []
                to_check = [ipix] ### add it to the list of things to check

                while len(to_check): # there are pixels in this mode we have to check
                        ipix = to_check.pop() # take one pixel from those to be checked.
			isinterior = True
                        for neighbour in hp.get_all_neighbours(nside, ipix, nest=nest):# get neighbors as rtheta, rphi

                                if neighbour == -1: ### when neighbour == -1, there is no corresponding pixel in this direction
                                        pass
				else:
					isinterior *= abstruth[neighbour] ### check to see if the neighbor is in the set
					                                  ### we don't care if it has already been visited.
	                                # try to find pixel in skymap
        	                        if truth[neighbour]: ### pixel in the set and has not been visited before
                	                        truth[neighbour] = False ### remove neighbour from global set
                                	        to_check.append( neighbour ) ### add to list of things to check
	                                else: ### pixel not in the set or has been visited before
        	                                pass
			if not isinterior:
				boarder.append( ipix )
                boarders.append( boarder )

        return boarders

###
def num_modes(posterior, theta, phi, nside=None, nest=False):
	"""
	computes the number of modes in the area bounded by theta, phi
	"""
	npix = len(posterior)
	if not nside:
		nside = hp.npix2nside(npix)
	pix = list(np.arange(npix)[posterior>=posterior[hp.ang2pix(nside, theta, phi, nest=nest)]]) ## get list of pixels
	return len( __into_modes(nside, pix, nest=nest) )

### 
def size_modes(posterior, theta, phi, nside=None, nest=False, degrees=False):
	"""
	computes the sizes of modes in the area bounded by theta, phi
	"""
	npix = len(posterior)
	if not nside:
		nside = hp.npix2nside(npix)
	pix = list(np.arange(npix)[posterior>=posterior[hp.ang2pix(nside, theta, phi, nest=nest)]])
	return [len(_)*hp.nside2pixarea(nside, degrees=degrees) for _ in __into_modes(nside, pix)]

def __into_modes(nside, pix, nest=False):
	"""
	divides the list of pixels (pix) into simply connected modes
	"""
	### establish an array representing the included pixels
	npix = hp.nside2npix(nside)

	truth = np.zeros((npix,),bool)
	truth[pix] = True

	pixnums = np.arange(npix) ### convenient array we establish once

	modes = []
	while truth.any():
		ipix = pixnums[truth][0] ### take the first pixel
		truth[ipix] = False ### remove it from the global set
		mode = [ipix]
		to_check = [ipix] ### add it to the list of things to check

		while len(to_check): # there are pixels in this mode we have to check
			ipix = to_check.pop() # take one pixel from those to be checked.

			for neighbour in hp.get_all_neighbours(nside, ipix, nest=nest):# get neighbors as rtheta, rphi

				if neighbour == -1: ### when neighbour == -1, there is no corresponding pixel in this direction
					pass
				# try to find pixel in skymap
				elif truth[neighbour]: ### pixel in the set and has not been visited before
					truth[neighbour] = False ### remove neighbour from global set
#					truth[neighbour] = 0 ### remove neighbour from global set
					mode.append( neighbour ) ### add to this mode
					to_check.append( neighbour ) ### add to list of things to check
				else: ### pixel not in the set or has been visited before
					pass 
		modes.append( mode )

	return modes

###
def entropy(posterior, base=2.0):
	"""
	computes the shannon entropy in the posterior
	we compute the entropy with base=base
		base=2 => bits
		base=e => nats
	"""
	posterior = posterior[posterior > 0.0] ### log(0)*0 -> 0
	return -np.sum( np.log(posterior)*posterior)/np.log(base)

###
def information(posterior, base=2.0):
	"""
	computes the shannon entropy in the posterior
	we compute the entropy with base=base
		base=2 => bits
		base=e => nats
	"""
	n = len(posterior)
	posterior = posterior[posterior > 0.0] ### log(0)*0 -> 0
	logbase = np.log(base)
	return (np.log(n) + np.sum( np.log(posterior)*posterior) )/np.log(base)

#=================================================
#
# methods for comparing two skymaps
#
#=================================================
def mse(posterior1, posterior2):
	"""
	computes the mean square error between the two posteriors
		sum (p1 - p2)**2
	"""
	return 1.0*np.sum( (posterior1-posterior2)**2 )/len(posterior1)

###
def peak_snr(posterior1, posterior2):
	"""
	computes the peak signal-to-noise ratio between the two posteriors
		max(posterior1)/mse , max(posterior2)/mse
	"""
	_mse = mse(posterior1, posterior2)**0.5
	return np.max(posterior1)/_mse, np.max(posterior2)/_mse

###
def fidelity(posterior1, posterior2):
	"""
	computes the fidelity between the two posteriors
		sum (p1*p2)**0.5
	"""
	return np.sum( (posterior1*posterior2)**0.5 )

###
def indep_joint_entropy(posterior1, posterior2, base=2.0):
	"""
	return entropy(posterior1*posterior2, base=base)
	"""
	joint = posterior1*posterior2
	joint /= np.sum(joint)
	return entropy(joint, base=base)

###
def indep_joint_entropy(posterior1, posterior2, base=2.0):
	"""
	return information(posterior1*posterior2, base=base)
	"""
	joint = posterior1*posterior2
	joint /= np.sum(joint)
	return information(joint, base=base)

###
def KLdivergence(posterior1, posterior2, base=2.0):
	"""
	computes the Kullback-Leibler divergence
		sum log(p1/p2) * p1
	"""
	truth = posterior1>0
	return np.sum( np.log(posterior1[truth]/posterior2[truth])*posterior1[truth] )/np.log(base)

###
def symmetric_KLdivergence(posterior1, posterior2, base=2.0):
	"""
	computes the symmetric Kullback-Leibler divergence
		sum log(p1/p2)*(p1 - p2)
	"""
	return KLdivergence(posterior1, posterior2, base=base) + KLdivergence(posterior2, posterior2, base=base)

###
def symmetric_KLdivergence_walk( posterior1, posterior2, base=2.0, nside=False, nest=False ):
	"""
	compute the symmetric Kullback-Leibler divergence, stepping down the resolution to address edge effects
		sum log(p1/p2)*(p1 - p2)
	returns the first finite symKL found
	"""
	if not nside:
		nside = hp.npix2nside( len(posterior1) )
	symKL = symmetric_KLdivergence( posterior1, posterior2, base=base )
	while (nside > 1) and (symKL == np.infty):
		nside = nside / 2
		posterior1 = resample(posterior1, nside, nest=nest)
		posterior2 = resample(posterior2, nside, nest=nest)
		symKL = symmetric_KLdivergence( posterior1, posterior2, base=base )

	return symKL, nside

###
def structural_similarity(posterior1, posterior2, c1=(0.01)**2, c2=(0.03)**2):
	"""
	computes the structural similarity
		(2*m1*m2+c1)*(2*v12 + c2) / (m1**2 + m2**2 + c1)*(v1 + v2 + c2)
	dynamic range of our pixels is 1
	m1 => mean(p1)
	v1 => var(p1)
	v12=> covar(p1,p2)
	"""
	m1 = np.mean(posterior1)
	m2 = np.mean(posterior2)
	covar = np.cov(posterior1, posterior2)
	return ( (2*m1*m2+c1)*(2*covar[0,1] + c2) ) / ( (m1**2 + m2**2 + c1)*(covar[0,0] + covar[1,1] + c2) )

### 
def pearson(posterior1, posterior2):
	"""
	computes the pearson correlation coefficient
	"""
	covar = np.cov(posterior1, posterior2)
	return covar[0,1]/(covar[0,0]*covar[1,1])**0.5

###
def dot(posterior1, posterior2):
	"""
	takes the "inner product" of the two posteriors and returns the cos(theta) between them
	"""
	return np.sum(posterior1*posterior2)/(np.sum(posterior1**2)*np.sum(posterior2**2))**0.5

###
def geometric_overlap(pix1, pix2, nside, degrees=False):
        """
        computes the amount of area in the intersection and union of confidence regions from p1 and p2 defined by p_value
        """
        npix = hp.nside2npix(nside)
        pixarea = hp.nside2pixarea(nside, degrees=degrees)

	posterior1 = np.zeros((npix,),int)
	posterior2 = np.zeros((npix,),int)

	posterior1[pix1] = 1
	posterior2[pix2] = 1

        intersection = np.sum( posterior1*posterior2 )
        return intersection*pixarea, (np.sum(posterior1+posterior2) - intersection)*pixarea

###
def spotcheck(posterior1, posterior2, conf):
	"""
	finds the confidence regions corresonding to conf in posterior1 and computes the contained probability in posterior2, and vice versa
	assumes posetrior1 and posterior2 have the same nside and ordering
        return np.sum( posterior2[ credible_region(posterior1, conf) ] ), np.sum( posterior1[ credible_region(posterior2, conf) ] )
	"""
	if isinstance(conf, (int,float)):
		conf = np.array([conf])
	return [np.sum(posterior2[pix]) for pix in credible_region(posterior1, conf) ], [np.sum(posterior1[pix]) for pix in credible_region(posterior2, conf) ]

###
def __fits_to_table(posterior):
	"""
	helper function that maps a FITS file into the smallest possible table
	for use in twoPt functions
	"""
	npix = len(posterior)
	nside = hp.npix2nside(npix)

        ind = np.arange(npix)[posterior>0] ### only take the pixels that will contribute to the correlation function
        t, p = hp.pix2ang(nside, ind)

	return np.transpose([t, p, posterior[ind]])

###
def twoPt_fitsfits(posterior1, posterior2, Nsamp=101):
	"""
	estimates the 2-pt correlation function between the 2 posteriors. 
	This function assumes both posteriors are HEALPix all-sky maps

	returns theta, 2ptCorr(theta) as np.arrays with theta in radians
	"""
	### extract parameters from the maps to set up iteration
	npix1 = len(posterior1)
	nside1 = hp.npix2nside(npix1)
	npix2 = len(posterior2)
	nside2 = hp.npix2nside(npix2)

	pixarea1 = hp.nside2pixarea(nside1)
	pixarea2 = hp.nside2pixarea(nside2)

	kde_bandwidth = (pixarea1+pixarea2) ### variance of KDE Gaussian

        ### delegate!
        return twoPt_tabletable(__fits_to_table(posterior1), __fits_to_table(posterior2), Nsamp=Nsamp, kde_bandwidth=kde_bandwidth)

###
def twoPt_fitsfits_fast(posterior1, posterior2, Nsamp=101):
        """
        estimates the 2-pt correlation function between the 2 posteriors. 
        This function assumes both posteriors are HEALPix all-sky maps

        returns theta, 2ptCorr(theta) as np.arrays with theta in radians

        WARNING: this builds a BIG matrix to try to compute this quickly and may be very memory intensive...
	    delegates to twoPt_tabletable_fast
        """
        ### extract parameters from the maps to set up iteration
        npix1 = len(posterior1)
        nside1 = hp.npix2nside(npix1)
        npix2 = len(posterior2)
        nside2 = hp.npix2nside(npix2)

        pixarea1 = hp.nside2pixarea(nside1)
        pixarea2 = hp.nside2pixarea(nside2)

        kde_bandwidth = (pixarea1+pixarea2)

	### delegate!
	return twoPt_tabletable_fast(__fits_to_table(posterior1), __fits_to_table(posterior2), Nsamp=Nsamp, kde_bandwidth=kde_bandwidth)

###
def twoPt_fitstable(posterior, table, Nsamp=101):
	"""
	estimates the 2-pt correlation function between 2 posteriors
	assumes "posterior" is a HEALPix all-sky map and "table" is a tabular data format
	table must be a numpy array of the form: [(theta, phi, weight), (theta, phi, weight), ...]
	"""
        ### extract parameters from the maps to set up iteration
        npix = len(posterior)
        nside = hp.npix2nside(npix)

        kde_bandwidth = hp.nside2pixarea(nside)

        ### delegate!
        return twoPt_tabletable(__fits_to_table(posterior), table, Nsamp=Nsamp, kde_bandwidth=kde_bandwidth)

###
def twoPt_fitstable_fast(posterior, table, Nsamp=101):
        """
        estimates the 2-pt correlation function between 2 posteriors
        assumes "posterior" is a HEALPix all-sky map and "table" is a tabular data format
        table must be a numpy array of the form: [(theta, phi, weight), (theta, phi, weight), ...]

        WARNING: this builds a BIG matrix to try to compute this quickly and may be very memory intensive...
        """
        ### extract parameters from the maps to set up iteration
        npix = len(posterior)
        nside = hp.npix2nside(npix)

        kde_bandwidth = hp.nside2pixarea(nside)

        ### delegate!
        return twoPt_tabletable_fast(__fits_to_table(posterior), table, Nsamp=Nsamp, kde_bandwidth=kde_bandwidth)

###
def twoPt_tabletable(table1, table2, kde_bandwidth=1.0, Nsamp=101):
	"""
	estimate the 2-pt correlation function between 2 posteriors
	assumes both are tabular data formats
	tables must be a numpy array of the form: [(theta, phi, weight), (theta, phi, weight), ...]

	kde_bandwidth is the vairance used in the Gaussian KDE
	"""
        N = (2*np.pi*kde_bandwidth)**-0.5
        n = 0.5/kde_bandwidth

        ### set up iteration
        theta = np.linspace(0, np.pi, Nsamp)
        count = np.zeros_like(theta, dtype=float)

        for t1, p1, post1 in table1:
                cosT1 = np.cos(t1)
                sinT1 = np.sin(t1)
                for t2, p2, post2 in table2:
                        dTheta = np.arccos(cosT1*np.cos(t2) + sinT1*np.sin(t2)*np.cos(p1-p2))
                        count += N*np.exp(-n*(dTheta-theta)**2) * post1*post2 ### normalize Gaussian by inner product

        return theta, count

###
def twoPt_tabletable_fast(table1, table2, kde_bandwidth=1.0, Nsamp=101):
        """
        estimate the 2-pt correlation function between 2 posteriors
        assumes both are tabular data formats
        tables must be a numpy array of the form: [(theta, phi, weight), (theta, phi, weight), ...]

	kde_bandwidth is the variance used in the Gaussian KDE

        WARNING: this builds a BIG matrix to try to compute this quickly and may be very memory intensive...
	"""
	### set up array indexing
	t1, p1, posterior1 = table1.transpose()
	t2, p2, posterior2 = table2.transpose()

	IND1, IND2 = np.meshgrid(np.arange(len(t1)), np.arange(len(t2)))
	IND1 = IND1.flatten()
	IND2 = IND2.flatten()

	### compute angular separation between all pairs of points
        dTheta = np.arccos( cos_dtheta(t1[IND1], p1[IND1], t2[IND2], p2[IND2]) )

	### compute sampling via KDE
	theta = np.linspace(0, np.pi, Nsamp)

        kernal = (2*np.pi*kde_bandwidth)**-0.5 * np.exp(-0.5*(np.outer(dTheta, np.ones_like(theta)) - np.outer(np.ones_like(dTheta), theta))**2/kde_bandwidth)
        weights = np.outer(posterior1[IND1]*posterior2[IND2], np.ones_like(theta))

        count = np.sum(kernal * weights, axis=0) ### sum over all pixel pairs

        return theta, count
