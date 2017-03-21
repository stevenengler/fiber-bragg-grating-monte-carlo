# Steven Engler
# Phys 3311, Lakehead University
# March 22, 2017
#
# This simulation code was written based on the paper:
#    Modeling Light Propagation Through Bragg Reflection Gratings With an Adapted Monte Carlo Method
#    Ben E. Pelleg and Adam K. Fontecchio
#    Journal of Lightwave Technology, Vol. 32, No. 12, June 15, 2014
#
import numpy as np
import pylab as plt
import random
#
def calc_A(n1, n2, w, c, d):
	speed_of_light = 3e8
	k1 = n1*w/speed_of_light
	k2 = n2*w/speed_of_light
	#
	first_term = np.cos(k2*d)
	second_term = (0.5)*1j*(k2/k1 + k1/k2)*np.sin(k2*d)
	factor = np.exp(1j*k1*c)
	#
	return factor*(first_term+second_term)
#
def calc_B(n1, n2, w, c, d):
	speed_of_light = 3e8
	k1 = n1*w/speed_of_light
	k2 = n2*w/speed_of_light
	#
	first_term = (0.5)*1j*(k2/k1 - k1/k2)*np.sin(k2*d)
	factor = np.exp(-1j*k1*c)
	#
	return factor*first_term
#
def calc_C(n1, n2, w, c, d):
	speed_of_light = 3e8
	k1 = n1*w/speed_of_light
	k2 = n2*w/speed_of_light
	#
	first_term = -(0.5)*1j*(k2/k1 - k1/k2)*np.sin(k2*d)
	factor = np.exp(1j*k1*c)
	#
	return factor*first_term
#
def calc_D(n1, n2, w, c, d):
	speed_of_light = 3e8
	k1 = n1*w/speed_of_light
	k2 = n2*w/speed_of_light
	#
	first_term = np.cos(k2*d)
	second_term = -(0.5)*1j*(k2/k1 + k1/k2)*np.sin(k2*d)
	factor = np.exp(-1j*k1*c)
	#
	return factor*(first_term+second_term)
#
def calc_prop_matrix(n1, n2, w, c, d):
	return np.array(
		[
		 [calc_A(n1, n2, w, c, d), calc_B(n1, n2, w, c, d)],
		 [calc_C(n1, n2, w, c, d), calc_D(n1, n2, w, c, d)]
		]
	)
#
def reflectance_s(theta_i, n1, n2):
	term_1 = n1*np.cos(theta_i)
	term_2 = n2*np.sqrt(1-(np.sin(theta_i)*(n1/n2))**2)
	#
	return ((term_1-term_2)/(term_1+term_2))**2
#
def reflectance_p(theta_i, n1, n2):
	term_1 = n1*np.sqrt(1-(np.sin(theta_i)*(n1/n2))**2)
	term_2 = n2*np.cos(theta_i)
	#
	return ((term_1-term_2)/(term_1+term_2))**2
#
def reflectance(theta_i, n1, n2):
	return 0.5*(reflectance_s(theta_i, n1, n2) + reflectance_p(theta_i, n1, n2))
#
def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
	return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)
#
if __name__ == '__main__':
	np.set_printoptions(precision=3)
	np.set_printoptions(suppress=True)
	#
	speed_of_light = 3e8
	#
	n1 = 1.5
	n2 = 1.7
	c = 75e-9
	d = c
	periods = 10 # this is half the num of periods they specify in the paper
	num_photons_to_simulate = 70
	#
	wavelengths = np.linspace(400e-9, 700e-9, 500)
	#
	reflected_photons = []
	transmitted_photons = []
	#
	reflected = []
	transmitted = []
	#
	for lmda_index in range(len(wavelengths)):
		lmda = wavelengths[lmda_index]
		w = 2*np.pi*speed_of_light/lmda
		#
		mat_n1_to_n2 = calc_prop_matrix(n1, n2, w, c, d)
		#
		E_field = np.zeros((periods*2,2), dtype=np.cfloat)
		E_field[-1] = [1, 0]
		#
		for x in range(len(E_field)-1,0,-1):
			E_field[x-1] = np.dot(mat_n1_to_n2, E_field[x])
		#
		num_of_photons = np.square(np.absolute(E_field))
		num_of_photons /= num_of_photons[0, 0]
		#
		reflected_photons.append(num_of_photons[0][1])
		transmitted_photons.append(num_of_photons[-1][0])
		#
		reflectivities = np.zeros((len(E_field)-1,2))
		for x in range(len(reflectivities)):
			if x == len(reflectivities)-1:
				Rbn = reflectance(0, n1, n2)
				Rfn = 1-Rbn
			else:
				Nfn = num_of_photons[x][0]
				Nbnp = num_of_photons[x+1][1]
				Nfnp = num_of_photons[x+1][0]
				#
				# if Rbn == 1
				Rfn_max = (Nfn + Nbnp*1 - Nfnp)/Nfn
				# if Rbn == 0
				Rfn_min = (Nfn + Nbnp*0 - Nfnp)/Nfn
				#
				# if Rfn == 1
				Rbn_max = (Nfnp + Nfn*1 - Nfn)/Nbnp
				# if Rfn == 0
				Rbn_min = (Nfnp + Nfn*0 - Nfn)/Nbnp
				#
				Rfn_max = min(Rfn_max, 1)
				Rfn_min = max(Rfn_min, 0)
				Rbn_max = min(Rbn_max, 1)
				Rbn_min = max(Rbn_min, 0)
				#
				Rfn = (Rfn_max+Rfn_min)/2
				Rbn = (Rbn_max+Rbn_min)/2
				#
				#print(Rfn, (Nfn + Nbnp*Rbn - Nfnp)/Nfn, isclose(Rfn, (Nfn + Nbnp*Rbn - Nfnp)/Nfn, abs_tol=0.001))
				#print(Rbn, (Nfnp + Nfn*Rfn - Nfn)/Nbnp, isclose(Rbn, (Nfnp + Nfn*Rfn - Nfn)/Nbnp, abs_tol=0.001))
				#
				assert isclose(Rfn, (Nfn + Nbnp*Rbn - Nfnp)/Nfn, abs_tol=0.001)
				assert isclose(Rbn, (Nfnp + Nfn*Rfn - Nfn)/Nbnp, abs_tol=0.001)
			#
			reflectivities[x] = [Rfn, Rbn]
		#
		num_reflected = 0
		num_transmitted = 0
		#
		for x in range(num_photons_to_simulate):
			state = [0, 0]
			# position 0, direction 0 (fwd)
			#
			while state != [0, 1] and state != [len(reflectivities), 0]:
				# while photon is not at start and not moving backwards
				#    and not at end and moving forwards
				#
				position = state[0]
				direction = state[1]
				#
				if direction == 0:
					# moving forwards
					reflect_probability = reflectivities[position,0]
					reflect = random.random()<=reflect_probability
				else:
					# moving backwards
					reflect_probability = reflectivities[position-1,1]
					reflect = random.random()<=reflect_probability
				#
				if reflect:
					direction = direction^1
				else:
					if direction == 0:
						# moving forwards
						position += 1
					else:
						# moving backwards
						position -= 1
					#
				#
				state = [position, direction]
			#
			#print('finished one', random.random())
			if state[1] == 1:
				num_reflected += 1
			else:
				num_transmitted += 1
			#
		#
		reflected.append(num_reflected/num_photons_to_simulate)
		transmitted.append(num_transmitted/num_photons_to_simulate)
	#
	plt.plot(wavelengths*1e9, reflected_photons, '.')
	plt.ylim([0, 1.05])
	plt.ylabel('Reflected Intensity')
	plt.xlabel('$\lambda$ (nm)')
	plt.show()
	#
	plt.plot(wavelengths*1e9, reflected, '.')
	plt.ylim([0, 1.05])
	plt.ylabel('Reflected Intensity')
	plt.xlabel('$\lambda$ (nm)')
	plt.show()
#
