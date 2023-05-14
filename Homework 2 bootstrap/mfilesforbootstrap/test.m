tt=[10,50,100,500];
fTstar=[0.0106, 0.0021, 0.0011, 8.2836e-05];

plot(tt,fTstar,'b','linewidth',2)
hold on
% legend('Non-param. bootstr.','Param. bootstrap','True distr. f_T(t)')
legend('exact bias from T_s_i_m_u_l_a_t_e')
set(gca,'fontsize',14,'fontweight','b')
xlabel('n_s_a_m_p_l_e')
ylabel('Bias')
print -dpdf 'exact bias from T_simulate.pdf'
hold off


