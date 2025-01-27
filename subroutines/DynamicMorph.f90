!	Name: DynamicMorph.f90

!	Purpose:
!    	This script is designed to compute displacements for the nodes within the ankle morph section of the ankle FE model using the UAMP subroutine.
!		Displacements are computed at the beginning of each increment using extracted rotations about the subtalar and talocrucal axes, and applied to
!		the nodes within the morph as an amplitude through "morph_include_20dfx_30sup_bcs.inc".

!	Usage:
!   	abaqus -j jobname user=DynamicMorph.f90
!		*Note* Above syntax must be submitted through the command line of the Intel Fortran Compiler.

!	Dependencies:
!   	- Requires "MorphData.inc" and "morph_include_20dfx_30sup_bcs.inc" to be placed in same parent directory.
!		- Requires access to Intel Fortran Compiler.

!	Authors:
!   	Noah Palermo (npalermo@mines.edu)
!		Mykola Mazur (mazur@mines.edu)

!	Last Modified:
!		August 25, 2024

!	Notes:
!		Ensure that this script is in the same parent directory as the CAE database.



!	Define global variables

module global_vars
	integer :: nodecounter = 1
	integer :: firstdof = 0
	integer :: dofcounter = 0
	integer :: inccounter = 1
	real :: currentnodedisp(3), a_shift(787, 3), a_shiftold(787, 3)
	logical :: firstinc = .true.
	logical :: newinc = .true.
	logical :: firstincrement = .true.
	common /counter_block/ counter
end module global_vars

!	** Begin subroutine UAMP **
SUBROUTINE UAMP(ampName, time, ampValueOld, dt, nProps, props, nSvars, &
          svars, lFlagsInfo, &
          nSensor, sensorValues, sensorNames, jSensorLookUpTable, &
          AmpValueNew, &
          lFlagsDefine, &
          AmpDerivative, AmpSecDerivative, AmpIncIntegral, &
          AmpDoubleIntegral, JInc)

	  USE global_vars
      INCLUDE 'ABA_PARAM.INC'


	!	Time indices
	parameter (iStepTime        = 1, &
                iTotalTime       = 2, &
                nTime            = 2)
	!   Flags passed in for information
	parameter (iInitialization   = 1, &
			iRegularInc       = 2, &
			iCuts             = 3, &
			ikStep            = 4, &
			nFlagsInfo        = 4)
	!   Optional flags to be defined
	parameter (iComputeDeriv       = 1, &
			iComputeSecDeriv    = 2, &
			iComputeInteg       = 3, &
			iComputeDoubleInteg = 4, &
			iStopAnalysis       = 5, &
			iConcludeStep       = 6, &
			nFlagsDefine        = 6)
	dimension time(nTime), lFlagsInfo(nFlagsInfo), &
		   lFlagsDefine(nFlagsDefine)
	dimension jSensorLookUpTable(*)
	dimension sensorValues(nSensor), svars(nSvars), props(nProps)
	character*80 sensorNames(nSensor)
	character*80 ampName
	real :: tc_rot, st_rot
	  

!	Declare morph variables

	!	Mixed array/real data structure type for morph output
	type :: MixedVector
		real, dimension(3, 3) :: matrix
		real :: real1
		real :: real2
		real :: real3
	end type MixedVector

	real :: testvector(3)
	real :: subdo(3), tc_ax(3), st_ax(3), tc_p1(3), tc_p2(3), JCSor(3), st_p1(3), st_p2(3)
	real :: p1(3), p2(3), p3(3), x(100), dfx(100), sup(100)
	real :: arot(787, 3), ashift(787, 3), amorph(787, 3), rad1(787, 1), a_rot(787, 3), rawankle(787, 3), proxarray(2, 3), distarray(2, 3)
	real :: Rst(3, 3), Rtc(3, 3)
	real :: r1(3), r2(3), r2_rot(3), c(3), cax(2, 3), r3(3), r3_rot(3), nca(3), pd(69), dd(92)
	real :: extrude, var, dfxi, supi, maxc, minc, popt_d, dopt_d, p_toa_dist, p_tod_dist
	real :: input1(3), input2(3), randvec1(3), randvec2(3), s
	integer :: i, j, n, k, popt_i, dopt_i, tempt_i(1), tempd_i(1), tempi_p(1), tempi_d(1)
	logical :: found_in_distal, found_in_proximal, rstzeros, rtczeros
	real :: a_morph(787, 3)
	type(MixedVector) :: output1, output2
	
	INCLUDE 'MorphData.inc'
	
	!	Get talocrural and subtalar axis rotation values for current increment (actually from past increment since UAMP is called
	!	at the beginning of each increment
	tc_rot = real(sensorValues(2))
	st_rot = real(sensorValues(1))
	
	!	Skip first time subroutine is called
	if (.not. firstinc) then
	
		!	Increment dof counter
		dofcounter = dofcounter + 1
		
		!	Check to see if job is on a new increment. If so, morph calculations will be done
		if (newinc) then
			if (.not. firstincrement) then
			!	Flip newinc off
			newinc = .false.
			
			!	** BEGIN MORPH SCRIPT **
			
			!	Define matrix comprised of all raw ankle data.
			do i = 1, 787
				rawankle(i, 1) = nodesankle(i, 2)
				rawankle(i, 2) = nodesankle(i, 3)
				rawankle(i, 3) = nodesankle(i, 4)
			end do

			!	Displacement from tib origin (JCS origin to subtalar origin)
			subdo(1) = -0.051042192975949 * 1000
			subdo(2) = -0.0439044493611044 * 1000
			subdo(3) = 0.00828899258498085 * 1000

			!	Direction of the talocrural axis defined in the tib reference frame
			tc_ax(1) = -0.10501355
			tc_ax(2) = -0.17402245
			tc_ax(3) = 0.97912632

			!	Direction of the subtalar axis defined in the tib reference frame
			st_ax(1) = 0.78717961
			st_ax(2) = 0.60474746
			st_ax(3) = -0.12094949

			!	Points defining the talocrural axis
			tc_p1(1) = -40.004421
			tc_p1(2) = 12.561185
			tc_p1(3) = 33.373463

			tc_p2 = tc_p1 - 75.25*tc_ax

			!	JCS origin (midpoint of tc_p1 and tc_p2)
			JCSor(1) = -36.053311
			JCSor(2) = 19.108792
			JCSor(3) = -3.466168

			!	Points defining the subtalar axis
			st_p1 = JCSor + subdo
			st_p2 = st_p1 + 90*st_ax

			!	Define talus
			p1 = JCSor
			p3 = st_p2
			p2(1) = p3(1)
			p2(2) = p1(2)
			p2(3) = p3(3)

			extrude = -(p3(3) - p1(3))

			!	Load node data for distal nodes
			do i = 1, 100
				x(i) = (i-1) / real(99)
			end do

			do i = 1, 100
				dfx(i) = cos(x(i)*2*3.14159265358979323846)
				sup(i) = sin(x(i)*2*3.14159265358979323846)
			end do

			!	Rotations pulled from Abaqus job
			dfxi = tc_rot
			supi = st_rot

			!	Define inputs for helical2cardan
			input1(1) = dfxi*tc_ax(1)
			input1(2) = dfxi*tc_ax(2)
			input1(3) = dfxi*tc_ax(3)

			input2(1) = supi*st_ax(1)
			input2(2) = supi*st_ax(2)
			input2(3) = supi*st_ax(3)

			!	Find rotation matrices for both inputs
			output1 = helical2cardan(input1)
			output2 = helical2cardan(input2)

			!	Define tc rotation matrix
			Rtc = output1%matrix


			!	Define st rotation matrix
			Rst = output2%matrix

			!	Fixed reference point on the talocrural axis
			r1 = tc_p1

			!	Rotating reference point on the subtalar
			r2 = st_p1 - r1
			r2_rot = matmul(Rtc, r2)

			!	Center point for arc that cuts proximal ankle surf
			c(1) = -34.78666
			c(2) = 53.87478
			c(3) = (maxval(rawankle(:, 3)) + minval(rawankle(:, 3))) / 2

			cax(1, 1) = -34.78666
			cax(1, 2) = 53.87478
			cax(1, 3) = minval(rawankle(:, 3))
			cax(2, 1) = -34.78666
			cax(2, 2) = 53.87478
			cax(2, 3) = maxval(rawankle(:, 3))

			!	Compute the full transformation due to dorsiflexion and supination for all nodes
			do i = 1, 787

				do n = 1, 3
					r3(n) = 0
					r3_rot(n) = 0
					nca(n) = 0
					a_rot(i, n) = 0
				end do

				r3 = rawankle(i, :) - st_p1
				r3_rot = matmul(matmul(Rtc, Rst), r3)

				a_rot(i, :) = r1 + r2_rot + r3_rot

				!	This is the unit vector pointing from the center to each vertex
				nca = findunit(rawankle(i, :) - c)

				!	Search through all proximal points and create an array containing their distance to the current
				!	vertex along the nca trajectory
				do n = 1, 69
					proxarray(1, :) = rawankle(i, :)
					proxarray(2, :) = nodesproximal(n, 2:4)
					pd(n) = pdist(proxarray)
				end do
				
				!	Search through all distal points and create an array containing their distance to the current
				!	vertex along the nca trajectory
				do n = 1, 92
					distarray(1, :) = rawankle(i, :)
					distarray(2, :) = nodesdistal(n, 2:4)
					dd(n) = pdist(distarray)
				end do
							 
				!	Which proximal point is closest to the current vertex?
				tempi_p(1) = 0
				popt_d = minval(pd)
				tempi_p = minloc(pd)
				popt_i = tempi_p(1)
				
				!	Which distal point is closest to the current vertex?
				tempi_d(1) = 0
				dopt_d = minval(dd)
				tempi_d = minloc(dd)
				dopt_i = tempi_d(1)
				
				!	What is the distance from proximal to distal along the proximal to distal path?
				p_toa_dist = dot(rawankle(i, :) - nodesproximal(popt_i, 2:4), findunit(nodesdistal(dopt_i, 2:4) - nodesproximal(popt_i, 2:4)))
				p_tod_dist = dot(nodesdistal(dopt_i, 2:4) - nodesproximal(popt_i, 2:4), findunit(nodesdistal(dopt_i, 2:4) - nodesproximal(popt_i, 2:4)))

				!	Compute scale factor s
				s = scales(p_toa_dist, p_tod_dist, 2)
				
				!	Check to see if current node is proximal or distal
				found_in_distal = any(nodesdistal(:, 1) == nodesankle(i, 1))
				found_in_proximal = any(nodesproximal(:, 1) == nodesankle(i, 1))
				
				if (found_in_distal) then
					!	Set s to 1 if found in distal
					s = 1.0
					
				else if (found_in_proximal) then
					!	Set s to 0 if found in proximal
					s = 0.0
					
				end if
				
				!	Perform shift calculation
				do k = 1, 3
					a_shift(i, k) = s*(a_rot(i, k) - rawankle(i, k))
					a_morph(i, k) = rawankle(i, k) + a_shift(i, k)
				end do
			
			end do
		
			else
				!	Set all shift values equal to zero
				do n = 1, 787
					do k = 1, 3
						a_shift(n, k) = 0
					end do
					
				end do
				
			end if
			
		end if
			
	
		

		!	** OPTIONAL PRINT STATEMENTS FOR DEBUGGING **
		!print*, '*********'
		!print*, 'current inc'
		!print*, inccounter
		!print*, 'tc rot'
		!print*, tc_rot
		!print*, 'st rot'
		!print*, st_rot
		!print*, 'current node'
		!print*, nodecounter
		!print*, 'current dof'
		!print*, dofcounter


		!	Define the output amplitude of UAMP as the shift for the current node and current dof
		AmpValueNew = a_shift(nodecounter, dofcounter)
		
		!	More debugging
		!print*, 'Current shift'
		!print*, AmpValueNew

		if (firstinc) then
			AmpValueNew = 0
		end if
		
		!	When the 3rd dof is reached, reset dof counter and increment node
		if (dofcounter == 3) then
			nodecounter = nodecounter + 1
			dofcounter = 0
		end if
		
		!	When the last node is reached, reset node counter and increment inc counter
		if (nodecounter == 788) then
			nodecounter = 1
			inccounter = inccounter + 1
			newinc = .true.
			firstincrement = .false.
		end if
	end if
	
	
	!	Flip firstinc after first time subroutine is called
	firstinc = .false.

	
    RETURN
	  
	  
	contains

		!	pdist resturns the linear distance between two points in cartesian space
		function pdist(mat) result(dist)
			real, intent(in) :: mat(2, 3)
			real :: dist
			dist = sqrt((mat(1, 1) - mat(2, 1)) ** 2 + (mat(1, 2) - mat(2, 2)) ** 2 + (mat(1, 3) - mat(2, 3)) ** 2)
		end function pdist

		!	norm computes the magnitude of an input vector
		function norm(vec) result(mag)
			real, intent(in) :: vec(3)
			real :: mag
			mag = sqrt(vec(1) ** 2 + vec(2) ** 2 + vec(3) ** 2)
		end function norm

		!	findunit finds the unit vector in the direction of an input vector
		function findunit(vec) result(unitvec)
			real, intent(in) :: vec(3)
			real :: unitvec(3)
			unitvec(1) = vec(1)/norm(vec)
			unitvec(2) = vec(2)/norm(vec)
			unitvec(3) = vec(3)/norm(vec)
		end function findunit
		
		!	dotprod finds the dot product of two input vectors
		function dot(vecA, vecB) result(dotprod)
			real, intent(in) :: vecA(3), vecB(3)
			real :: dotprod
			dotprod = vecA(1)*vecB(1) + vecA(2)*vecB(2) + vecA(3)*vecB(3)
		end function dot

		!	cross finds the cross
		function cross(vecA, vecB) result(crossprod)
		 real, intent(in) :: vecA(3), vecB(3)
		 real :: crossprod(3)
		 crossprod(1) = vecA(2) * vecB(3) - vecA(3) * vecB(2)
		 crossprod(2) = vecA(3) * vecB(1) - vecA(1) * vecB(3)
		 crossprod(3) = vecA(1) * vecB(2) - vecA(2) * vecB(1)
		end function cross
		
		function findtranspose(matA) result(transmat)
		 real, intent(in) :: matA(3, 3)
		 real :: transmat(3, 3)
		 integer :: i, j

		 ! Transpose the matrix
		 do i = 1, 3
		  do j = 1, 3
		   transmat(j, i) = matA(i, j)
		  end do
		 end do

		end function findtranspose
		
	!	function findvectranspose(vec) result(colvec)
	!	 real, intent(in) :: vec(3)
	!	 real :: colvec(1, 3)
		
		function matbyvec(mat, vec) result(outvec)
		 real, intent(in) :: mat(3, 3), vec(3)
		 real :: outvec(3)
		 integer :: i, j
		 
		 do i = 1, 3
		  do j = 1, 3
		   outvec(i) = outvec(i) + mat(i, j) * vec(i)
		  end do
		 end do
		end function matbyvec
		 
		function matrixmult(mat, vec) result(outvec)
		 real, intent(in) :: mat(3, 3), vec(3)
		 real :: outvec(3)
		 integer :: i, j
		 
		 do i = 1, 3
		  do j = 1, 3
		   outvec(i) = outvec(i) + mat(i, j) * vec(j)
		  end do
		 end do
		end function matrixmult
		
		function elementwisemult(matA, matB) result(outmat)
		 real, intent(in) :: matA(3, 3), matB(3, 3)
		 real :: outmat(3, 3)
		 integer :: i, j, k
		 do i = 1, 3
		  do j = 1, 3
		   outmat(i, j) = 0
		  end do
		 end do
		 
		 do i = 1, 3
		  do j = 1, 3
		   do k = 1, 3
			outmat(i, j) = outmat(i, j) + matA(i, k) * matB(k, j)
		   end do
		  end do
		 end do
		 
		end function elementwisemult
		
		function scales(in1, in2, method) result(s)
		 real, intent(in) :: in1, in2
		 integer, intent(in) :: method
		 real :: s
		 
		 if (method == 1) then
		  s = in1/in2
		 else if (method == 2) then
		  s = sin(3.14159265358979323846*(in1/in2)/2)
		 end if
		end function scales
		 

	function helical2cardan(alpha_vector) result(card)
	 real, intent(in) :: alpha_vector(3)
	 real :: alpha, i, flx, lat, axl, j
	 real :: p(3), hx(3), hy(3), hz(3), yvec(3), zvec(3), L(3), randvec(3)
	 real :: Rx(3, 3), Rxt(3, 3), xlocal(3), xrot(3), hxyz(3, 3), hxyzt(3, 3), bx(3), R(3, 3)
	 real :: ylocal(3), yrot(3), by(3), zlocal(3), zrot(3), bz(3)
	 real :: Rx1(3, 3), Ry1(3, 3), Rz1(3, 3), Ri(3, 3)
	 
	 type(MixedVector) :: card
	 

	 
	 ! Initialize variables to ensure no rollover from previous use of function
	 
	 do i = 1, 3
	  p(i) = 0
	  hx(i) = 0
	  hy(i) = 0
	  hz(i) = 0
	  yvec(i) = 0
	  zvec(i) = 0
	  L(i) = 0
	  xlocal(i) = 0
	  xrot(i) = 0
	  bx(i) = 0
	  ylocal(i) = 0
	  yrot(i) = 0
	  by(i) = 0
	  zlocal(i) = 0
	  zrot(i) = 0
	  bz(i) = 0
	  randvec(i) = 0
	 end do
	 
	 alpha = 0
	 flx = 0
	 lat = 0
	 axl = 0
	 
	 
	 do i = 1, 3
	  do j = 1, 3
	   Rx(i, j) = 0
	   Rxt(i, j) = 0
	   hxyz(i, j) = 0
	   hxyzt(i, j) = 0
	   R(i, j) = 0
	   Rx1(i, j) = 0
	   Ry1(i, j) = 0
	   Rz1(i, j) = 0
	   Ri(i, j) = 0
	  end do
	 end do
	 
	 
	 ! Compute magnitude of input vector
	 alpha = norm(alpha_vector)
	 
	 ! Rotation vector components are used to define the rotation axis
	 p(1) = alpha_vector(1)/alpha
	 p(2) = alpha_vector(2)/alpha
	 p(3) = alpha_vector(3)/alpha

	 
	 ! Make arbitrary local frame with x-axis lying along helical rotation axis
	 hx = p
	 
	 ! Create random vector
		 
	 ! Initialize random seed
	 
	 call random_seed()

	 call random_number(randvec)
		 
	 ! Find part of random vector perpendicular to p (unit vector)
	 hy = findunit(randvec - dot(randvec, p)*p)
	 
	 ! Simple cross product for z-axis definition
	 hz = cross(hx, hy)
	 

	 ! Define matrix composed of all three axes
	 hxyz(1, 1) = hx(1)
	 hxyz(1, 2) = hx(2)
	 hxyz(1, 3) = hx(3)
	 hxyz(2, 1) = hy(1)
	 hxyz(2, 2) = hy(2)
	 hxyz(2, 3) = hy(3)
	 hxyz(3, 1) = hz(1)
	 hxyz(3, 2) = hz(2)
	 hxyz(3, 3) = hz(3)

	 
	 ! Find local representation of vector we wish to rotate by projecting it into
	 ! the arbitrary local frame we constructed above. We do this for x, y, and z
	 
	 ! X:
	 
	 xlocal(1) = hx(1)
	 xlocal(2) = hy(1)
	 xlocal(3) = hz(1)
	 
	 ! Define rotation matrix Rx
	 Rx(1, 1) = 1
	 Rx(1, 2) = 0
	 Rx(1, 3) = 0
	 Rx(2, 1) = 0
	 Rx(2, 2) = cos(alpha)
	 Rx(2, 3) = sin(alpha)
	 Rx(3, 1) = 0
	 Rx(3, 2) = -sin(alpha)
	 Rx(3, 3) = cos(alpha)
	 
	 ! Rotate xlocal around the rotation axis, which is local x axis; note that
	 ! we use the transpose of Rx since we are rotating the vector rather than the
	 ! reference frame
	 Rxt = findtranspose(Rx)
	 
	 xrot = matmul(Rxt, xlocal)
	 
	 ! Find transpose of our hx, hy, hz matrix
	 hxyzt = findtranspose(hxyz)
	 
	 ! Put xrot back into global components by reversing the local transformation
	 bx = matmul(hxyzt, xrot)
	 
	 
	 ! Y:
	 
	 ylocal(1) = hx(2)
	 ylocal(2) = hy(2)
	 ylocal(3) = hz(2)
	 
	 ! Rotate ylocal around the rotation axis, which is the local x axis; note that
	 ! we use the transpose of Rx since we are rotating the vector rather than the
	 ! reference frame
	 yrot = matmul(Rxt, ylocal)
	 
	 
	 ! Put yrot back into global components by reversing the local transformation
	 by = matmul(hxyzt, yrot)
	 

	 ! Z:
	 
	 zlocal(1) = hx(3)
	 zlocal(2) = hy(3)
	 zlocal(3) = hz(3)
	 
	 ! Rotate zlocal around the rotation axis, which is the local x axis; note that
	 ! we use the transpose of Rx since we are rotating the vector rather than the
	 ! reference frame	 
	 zrot = matmul(Rxt, zlocal)
	 
	 ! Put zrot back into global components by reversing the local transformation
	 bz = matmul(hxyzt, zrot)

	 ! Define Cardan angles
	 
	 yvec(1) = 0
	 yvec(2) = 1
	 yvec(3) = 0
	 
	 zvec(1) = 0
	 zvec(2) = 0
	 zvec(3) = 1
	 
	 L = cross(by, zvec)
	 L = findunit(L)
	 flx = dot(L, yvec)
	 flx = asin(flx)
	 lat = dot(zvec, by)
	 lat = asin(lat)
	 axl = dot(L, bz)
	 axl = asin(axl)
	 
	 Rz1(1, 1) = cos(flx)
	 Rz1(1, 2) = sin(flx)
	 Rz1(1, 3) = 0
	 Rz1(2, 1) = -sin(flx)
	 Rz1(2, 2) = cos(flx)
	 Rz1(2, 3) = 0
	 Rz1(3, 1) = 0
	 Rz1(3, 2) = 0
	 Rz1(3, 3) = 1
	 
	 Rx1(1, 1) = 1
	 Rx1(1, 2) = 0
	 Rx1(1, 3) = 0
	 Rx1(2, 1) = 0
	 Rx1(2, 2) = cos(lat)
	 Rx1(2, 3) = sin(lat)
	 Rx1(3, 1) = 0
	 Rx1(3, 2) = -sin(lat)
	 Rx1(3, 3) = cos(lat)
	 
	 Ry1(1, 1) = cos(axl)
	 Ry1(1, 2) = 0
	 Ry1(1, 3) = -sin(axl)
	 Ry1(2, 1) = 0
	 Ry1(2, 2) = 1
	 Ry1(2, 3) = 0
	 Ry1(3, 1) = sin(axl)
	 Ry1(3, 2) = 0
	 Ry1(3, 3) = cos(axl)
	 
	 Ri = matmul(Ry1, Rx1)
	 R = matmul(Ri, Rz1)
	 
	 R = findtranspose(R)
	 
	 card%matrix = R
	 card%real1 = flx
	 card%real2 = lat
	 card%real3 = axl
	 
	 
	 
	end function helical2cardan
    END