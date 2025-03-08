	module constants
  		implicit none
  		integer, parameter :: ip = 4                                   ! 4-byte integer
		integer, parameter :: rp = selected_real_kind(p=15, r=307)     ! 8-byte (double-precision) real
	end module constants
