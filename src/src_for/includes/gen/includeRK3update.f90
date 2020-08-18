q (i,j,1) = param_float(fadrTIMESTEP)*rk1(nrk)*rhs(i,j,1) + q2(i,j,1)

q (i,j,2) = param_float(fadrTIMESTEP)*rk1(nrk)*rhs(i,j,2) + q2(i,j,2)

q (i,j,3) = param_float(fadrTIMESTEP)*rk1(nrk)*rhs(i,j,3) + q2(i,j,3)

q (i,j,4) = param_float(fadrTIMESTEP)*rk1(nrk)*rhs(i,j,4) + q2(i,j,4)

one_over_rho = 1.0_wp/q(i,j,1)
q(i,j,2) = q(i,j,2)*one_over_rho
q(i,j,3) = q(i,j,3)*one_over_rho
q(i,j,4) = q(i,j,4)*one_over_rho
