function [Voltage_output, phasors] = ieee34_iter2(loads, Zbus, Cbus)

tapA814R = 1.0 + 0.00625 * 12;
tapB814R = 1.0 + 0.00625 * 5;
tapC814R = 1.0 + 0.00625 * 5;
alpha814R = [tapA814R; tapB814R; tapC814R];
alphaM814R = alpha814R * alpha814R';

tapA852R = 1.0 + 0.00625 * 13;
tapB852R = 1.0 + 0.00625 * 11;
tapC852R = 1.0 + 0.00625 * 12;
alpha852R = [tapA852R; tapB852R; tapC852R];
alphaM852R = alpha852R * alpha852R';

Xfmr_ratio = 4.16/24.9;
Xfmr = [4.16/24.9; 4.16/24.9; 4.16/24.9];
XfmrM = Xfmr * Xfmr';

Z800802 = Zbus([4, 5, 6],[7, 8, 9]);
Z802806 = Zbus([7, 8, 9],[10, 11, 12]);
Z806808 = Zbus([10, 11, 12],[13, 14, 15]);
Z808810 = Zbus([14],[16]);
Z808812 = Zbus([13, 14, 15],[17, 18, 19]);
Z812814 = Zbus([17, 18, 19],[20, 21, 22]);
Z814R850 = Zbus([23, 24, 25],[26, 27, 28]);
Z816818 = Zbus([29],[30]);
Z816824 = Zbus([29, 31, 32],[33, 34, 35]);
Z818820 = Zbus([30],[36]);
Z820822 = Zbus([36],[37]);
Z824826 = Zbus([34],[38]);
Z824828 = Zbus([33, 34, 35],[39, 40, 41]);
Z828830 = Zbus([39, 40, 41],[42, 43, 44]);
Z830854 = Zbus([42, 43, 44],[45, 46, 47]);
Z832858 = Zbus([48, 49, 50],[51, 52, 53]);
Z832888 = Zbus([48, 49, 50],[84, 85, 86]);
Z834860 = Zbus([54, 55, 56],[57, 58, 59]);
Z834842 = Zbus([54, 55, 56],[60, 61, 62]);
Z836840 = Zbus([63, 64, 65],[66, 67, 68]);
Z836862 = Zbus([63, 64, 65],[69, 70, 71]);
Z842844 = Zbus([60, 61, 62],[72, 73, 74]);
Z844846 = Zbus([72, 73, 74],[75, 76, 77]);
Z846848 = Zbus([75, 76, 77],[78, 79, 80]);
Z850816 = Zbus([26, 27, 28],[29, 31, 32]);
Z852R832 = Zbus([81, 82, 83],[48, 49, 50]);
Z854856 = Zbus([46],[87]);
Z854852 = Zbus([45, 46, 47],[88, 89, 90]);
Z858864 = Zbus([51],[91]);
Z858834 = Zbus([51, 52, 53],[54, 55, 56]);
Z860836 = Zbus([57, 58, 59],[63, 64, 65]);
Z862838 = Zbus([70],[92]);
Z888890 = Zbus([84, 85, 86],[93, 94, 95]);
ZSOURCEBUS800 = Zbus([1, 2, 3],[4, 5, 6]);
Z814814R = Zbus([20, 21, 22],[23, 24, 25]);
Z852852R = Zbus([88, 89, 90],[81, 82, 83]);


% three phase voltage at slack bus
Vbase = 24900 / sqrt(3);
v0=1.00 * Vbase * [0,sqrt(3),0]';
% voltage upper and lower bounds
V_lb = 0.80 * Vbase;
V_ub = 1.10 * Vbase;
v_lb = V_lb * V_lb;
v_ub = V_ub * V_ub;

% sequential component parameters
a = -0.5 + 0.5 * i * sqrt(3);
A = 1/sqrt(3) * [1,1,1; 1, a*a, a; 1, a, a*a];
AH = 1/sqrt(3) * [1,1,1; 1, a, a*a; 1, a*a, a];



cvx_begin sdp quiet
% cvx_begin sdp 
% the solver: 
% cvx_solver SeDuMi;
cvx_solver Mosek;

% voltage square variables
variable v802(3,3) hermitian
variable v806(3,3) hermitian
variable v808(3,3) hermitian
variable v810(1,1) hermitian
variable v808_abc(3,3) hermitian
variable v812(3,3) hermitian
variable v814(3,3) hermitian
variable v850(3,3) hermitian
variable v818(1,1) hermitian
variable v816_abc(3,3) hermitian
variable v824(3,3) hermitian
variable v820(1,1) hermitian
variable v822(1,1) hermitian
variable v826(1,1) hermitian
variable v824_abc(3,3) hermitian
variable v828(3,3) hermitian
variable v830(3,3) hermitian
variable v854(3,3) hermitian
variable v858(3,3) hermitian
variable v888(3,3) hermitian
variable v860(3,3) hermitian
variable v842(3,3) hermitian
variable v840(3,3) hermitian
variable v862(3,3) hermitian
variable v844(3,3) hermitian
variable v846(3,3) hermitian
variable v848(3,3) hermitian
variable v816(3,3) hermitian
variable v832(3,3) hermitian
variable v856(1,1) hermitian
variable v854_abc(3,3) hermitian
variable v852(3,3) hermitian
variable v864(1,1) hermitian
variable v858_abc(3,3) hermitian
variable v834(3,3) hermitian
variable v836(3,3) hermitian
variable v838(1,1) hermitian
variable v862_abc(3,3) hermitian
variable v890(3,3) hermitian
variable v800(3,3) hermitian
variable v814R(3,3) hermitian
variable v852R(3,3) hermitian
variable vSOURCEBUS(3,3) hermitian

% complex power variables
variable S800802(3,3) complex
variable S802806(3,3) complex
variable S806808(3,3) complex
variable S808810(1,1) complex
variable S808812(3,3) complex
variable S812814(3,3) complex
variable S814R850(3,3) complex
variable S816818(1,1) complex
variable S816824(3,3) complex
variable S818820(1,1) complex
variable S820822(1,1) complex
variable S824826(1,1) complex
variable S824828(3,3) complex
variable S828830(3,3) complex
variable S830854(3,3) complex
variable S832858(3,3) complex
variable S832888(3,3) complex
variable S834860(3,3) complex
variable S834842(3,3) complex
variable S836840(3,3) complex
variable S836862(3,3) complex
variable S842844(3,3) complex
variable S844846(3,3) complex
variable S846848(3,3) complex
variable S850816(3,3) complex
variable S852R832(3,3) complex
variable S854856(1,1) complex
variable S854852(3,3) complex
variable S858864(1,1) complex
variable S858834(3,3) complex
variable S860836(3,3) complex
variable S862838(1,1) complex
variable S888890(3,3) complex
variable SSOURCEBUS800(3,3) complex
variable S814814R(3,3) complex
variable S852852R(3,3) complex

% current square variables
variable l800802(3,3) hermitian
variable l802806(3,3) hermitian
variable l806808(3,3) hermitian
variable l808810(1,1) hermitian
variable l808812(3,3) hermitian
variable l812814(3,3) hermitian
variable l814R850(3,3) hermitian
variable l816818(1,1) hermitian
variable l816824(3,3) hermitian
variable l818820(1,1) hermitian
variable l820822(1,1) hermitian
variable l824826(1,1) hermitian
variable l824828(3,3) hermitian
variable l828830(3,3) hermitian
variable l830854(3,3) hermitian
variable l832858(3,3) hermitian
variable l832888(3,3) hermitian
variable l834860(3,3) hermitian
variable l834842(3,3) hermitian
variable l836840(3,3) hermitian
variable l836862(3,3) hermitian
variable l842844(3,3) hermitian
variable l844846(3,3) hermitian
variable l846848(3,3) hermitian
variable l850816(3,3) hermitian
variable l852R832(3,3) hermitian
variable l854856(1,1) hermitian
variable l854852(3,3) hermitian
variable l858864(1,1) hermitian
variable l858834(3,3) hermitian
variable l860836(3,3) hermitian
variable l862838(1,1) hermitian
variable l888890(3,3) hermitian
variable lSOURCEBUS800(3,3) hermitian
variable l814814R(3,3) hermitian
variable l852852R(3,3) hermitian


minimize(trace(real(A * Z800802*l800802 * AH)) + trace(imag(A * Z800802*l800802 * AH)) + trace(real(A * Z802806*l802806 * AH)) + trace(imag(A * Z802806*l802806 * AH)) + trace(real(A * Z806808*l806808 * AH)) + trace(imag(A * Z806808*l806808 * AH)) + trace(real(Z808810*l808810)) + trace(imag(Z808810*l808810)) + trace(real(A * Z808812*l808812 * AH)) + trace(imag(A * Z808812*l808812 * AH)) + trace(real(A * Z812814*l812814 * AH)) + trace(imag(A * Z812814*l812814 * AH)) + trace(real(A * Z814R850*l814R850 * AH)) + trace(imag(A * Z814R850*l814R850 * AH)) + trace(real(Z816818*l816818)) + trace(imag(Z816818*l816818)) + trace(real(A * Z816824*l816824 * AH)) + trace(imag(A * Z816824*l816824 * AH)) + trace(real(Z818820*l818820)) + trace(imag(Z818820*l818820)) + trace(real(Z820822*l820822)) + trace(imag(Z820822*l820822)) + trace(real(Z824826*l824826)) + trace(imag(Z824826*l824826)) + trace(real(A * Z824828*l824828 * AH)) + trace(imag(A * Z824828*l824828 * AH)) + trace(real(A * Z828830*l828830 * AH)) + trace(imag(A * Z828830*l828830 * AH)) + trace(real(A * Z830854*l830854 * AH)) + trace(imag(A * Z830854*l830854 * AH)) + trace(real(A * Z832858*l832858 * AH)) + trace(imag(A * Z832858*l832858 * AH)) + trace(real(A * Z832888*l832888 * AH)) + trace(imag(A * Z832888*l832888 * AH)) + trace(real(A * Z834860*l834860 * AH)) + trace(imag(A * Z834860*l834860 * AH)) + trace(real(A * Z834842*l834842 * AH)) + trace(imag(A * Z834842*l834842 * AH)) + trace(real(A * Z836840*l836840 * AH)) + trace(imag(A * Z836840*l836840 * AH)) + trace(real(A * Z836862*l836862 * AH)) + trace(imag(A * Z836862*l836862 * AH)) + trace(real(A * Z842844*l842844 * AH)) + trace(imag(A * Z842844*l842844 * AH)) + trace(real(A * Z844846*l844846 * AH)) + trace(imag(A * Z844846*l844846 * AH)) + trace(real(A * Z846848*l846848 * AH)) + trace(imag(A * Z846848*l846848 * AH)) + trace(real(A * Z850816*l850816 * AH)) + trace(imag(A * Z850816*l850816 * AH)) + trace(real(A * Z852R832*l852R832 * AH)) + trace(imag(A * Z852R832*l852R832 * AH)) + trace(real(Z854856*l854856)) + trace(imag(Z854856*l854856)) + trace(real(A * Z854852*l854852 * AH)) + trace(imag(A * Z854852*l854852 * AH)) + trace(real(Z858864*l858864)) + trace(imag(Z858864*l858864)) + trace(real(A * Z858834*l858834 * AH)) + trace(imag(A * Z858834*l858834 * AH)) + trace(real(A * Z860836*l860836 * AH)) + trace(imag(A * Z860836*l860836 * AH)) + trace(real(Z862838*l862838)) + trace(imag(Z862838*l862838)) + trace(real(A * Z888890*l888890 * AH)) + trace(imag(A * Z888890*l888890 * AH)) + trace(real(A * ZSOURCEBUS800*lSOURCEBUS800 * AH)) + trace(imag(A * ZSOURCEBUS800*lSOURCEBUS800 * AH)) + trace(real(A * Z814814R*l814814R * AH)) + trace(imag(A * Z814814R*l814814R * AH)) + trace(real(A * Z852852R*l852852R * AH)) + trace(imag(A * Z852852R*l852852R * AH)) + 0)
subject to


% constraints: 
% (1): voltage lower and upper bounds 
v_lb <= diag(A * v802 * ctranspose(A)) <= v_ub;
v_lb <= diag(A * v806 * ctranspose(A)) <= v_ub;
v_lb <= diag(A * v808 * ctranspose(A)) <= v_ub;
v_lb <= diag(v810) <= v_ub;
v_lb <= diag(A * v812 * ctranspose(A)) <= v_ub;
v_lb <= diag(A * v814 * ctranspose(A)) <= v_ub;
v_lb <= diag(A * v850 * ctranspose(A)) <= v_ub;
v_lb <= diag(v818) <= v_ub;
v_lb <= diag(A * v824 * ctranspose(A)) <= v_ub;
v_lb <= diag(v820) <= v_ub;
v_lb <= diag(v822) <= v_ub;
v_lb <= diag(v826) <= v_ub;
v_lb <= diag(A * v828 * ctranspose(A)) <= v_ub;
v_lb <= diag(A * v830 * ctranspose(A)) <= v_ub;
v_lb <= diag(A * v854 * ctranspose(A)) <= v_ub;
v_lb <= diag(A * v858 * ctranspose(A)) <= v_ub;
Xfmr_ratio^2 * v_lb <= diag(A * v888 * ctranspose(A)) <= Xfmr_ratio^2 * v_ub;
v_lb <= diag(A * v860 * ctranspose(A)) <= v_ub;
v_lb <= diag(A * v842 * ctranspose(A)) <= v_ub;
v_lb <= diag(A * v840 * ctranspose(A)) <= v_ub;
v_lb <= diag(A * v862 * ctranspose(A)) <= v_ub;
v_lb <= diag(A * v844 * ctranspose(A)) <= v_ub;
v_lb <= diag(A * v846 * ctranspose(A)) <= v_ub;
v_lb <= diag(A * v848 * ctranspose(A)) <= v_ub;
v_lb <= diag(A * v816 * ctranspose(A)) <= v_ub;
v_lb <= diag(A * v832 * ctranspose(A)) <= v_ub;
v_lb <= diag(v856) <= v_ub;
v_lb <= diag(A * v852 * ctranspose(A)) <= v_ub;
v_lb <= diag(v864) <= v_ub;
v_lb <= diag(A * v834 * ctranspose(A)) <= v_ub;
v_lb <= diag(A * v836 * ctranspose(A)) <= v_ub;
v_lb <= diag(v838) <= v_ub;
Xfmr_ratio^2 * v_lb <= diag(A * v890 * ctranspose(A)) <= Xfmr_ratio^2 * v_ub;
v_lb <= diag(A * v800 * ctranspose(A)) <= v_ub;
v_lb <= diag(A * v814R * ctranspose(A)) <= v_ub;
v_lb <= diag(A * v852R * ctranspose(A)) <= v_ub;
vSOURCEBUS == v0 * ctranspose(v0);

% (1): voltage across a line 
v802 == v800 - S800802*ctranspose(Z800802) - Z800802*ctranspose(S800802) + Z800802*l800802*ctranspose(Z800802);
[v800, S800802; ctranspose(S800802), l800802] >= 0;
diag(A *(S800802-Z800802*l800802) * AH)- loads([7, 8, 9]) + diag(A * v802 * Cbus([7, 8, 9],[7, 8, 9]) * AH) == diag(A * S802806 * AH) + 0;

v806 == v802 - S802806*ctranspose(Z802806) - Z802806*ctranspose(S802806) + Z802806*l802806*ctranspose(Z802806);
[v802, S802806; ctranspose(S802806), l802806] >= 0;
diag(A *(S802806-Z802806*l802806) * AH)- loads([10, 11, 12]) + diag(A * v806 * Cbus([10, 11, 12],[10, 11, 12]) * AH) == diag(A * S806808 * AH) + 0;

v808 == v806 - S806808*ctranspose(Z806808) - Z806808*ctranspose(S806808) + Z806808*l806808*ctranspose(Z806808);
[v806, S806808; ctranspose(S806808), l806808] >= 0;
diag(A *(S806808-Z806808*l806808) * AH)- loads([13, 14, 15]) + diag(A * v808 * Cbus([13, 14, 15],[13, 14, 15]) * AH) == [0; diag(S808810); 0] + diag(A * S808812 * AH) + 0;

v810 == v808_abc([2],[2]) - S808810*ctranspose(Z808810) - Z808810*ctranspose(S808810) + Z808810*l808810*ctranspose(Z808810);
[v808_abc([2],[2]), S808810; ctranspose(S808810), l808810] >= 0;
diag(S808810-Z808810*l808810)- loads([16]) + diag(v810 * Cbus([16],[16])) == 0;

v812 == v808 - S808812*ctranspose(Z808812) - Z808812*ctranspose(S808812) + Z808812*l808812*ctranspose(Z808812);
[v808, S808812; ctranspose(S808812), l808812] >= 0;
diag(A *(S808812-Z808812*l808812) * AH)- loads([17, 18, 19]) + diag(A * v812 * Cbus([17, 18, 19],[17, 18, 19]) * AH) == diag(A * S812814 * AH) + 0;

v814 == v812 - S812814*ctranspose(Z812814) - Z812814*ctranspose(S812814) + Z812814*l812814*ctranspose(Z812814);
[v812, S812814; ctranspose(S812814), l812814] >= 0;
diag(A *(S812814-Z812814*l812814) * AH)- loads([20, 21, 22]) + diag(A * v814 * Cbus([20, 21, 22],[20, 21, 22]) * AH) == diag(A * S814814R * AH) + 0;

v850 == v814R - S814R850*ctranspose(Z814R850) - Z814R850*ctranspose(S814R850) + Z814R850*l814R850*ctranspose(Z814R850);
[v814R, S814R850; ctranspose(S814R850), l814R850] >= 0;
diag(A *(S814R850-Z814R850*l814R850) * AH)- loads([26, 27, 28]) + diag(A * v850 * Cbus([26, 27, 28],[26, 27, 28]) * AH) == diag(A * S850816 * AH) + 0;

v818 == v816_abc([1],[1]) - S816818*ctranspose(Z816818) - Z816818*ctranspose(S816818) + Z816818*l816818*ctranspose(Z816818);
[v816_abc([1],[1]), S816818; ctranspose(S816818), l816818] >= 0;
diag(S816818-Z816818*l816818)- loads([30]) + diag(v818 * Cbus([30],[30])) == diag(S818820) + 0;

v824 == v816 - S816824*ctranspose(Z816824) - Z816824*ctranspose(S816824) + Z816824*l816824*ctranspose(Z816824);
[v816, S816824; ctranspose(S816824), l816824] >= 0;
diag(A *(S816824-Z816824*l816824) * AH)- loads([33, 34, 35]) + diag(A * v824 * Cbus([33, 34, 35],[33, 34, 35]) * AH) == [0; diag(S824826); 0] + diag(A * S824828 * AH) + 0;

v820 == v818 - S818820*ctranspose(Z818820) - Z818820*ctranspose(S818820) + Z818820*l818820*ctranspose(Z818820);
[v818, S818820; ctranspose(S818820), l818820] >= 0;
diag(S818820-Z818820*l818820)- loads([36]) + diag(v820 * Cbus([36],[36])) == diag(S820822) + 0;

v822 == v820 - S820822*ctranspose(Z820822) - Z820822*ctranspose(S820822) + Z820822*l820822*ctranspose(Z820822);
[v820, S820822; ctranspose(S820822), l820822] >= 0;
diag(S820822-Z820822*l820822)- loads([37]) + diag(v822 * Cbus([37],[37])) == 0;

v826 == v824_abc([2],[2]) - S824826*ctranspose(Z824826) - Z824826*ctranspose(S824826) + Z824826*l824826*ctranspose(Z824826);
[v824_abc([2],[2]), S824826; ctranspose(S824826), l824826] >= 0;
diag(S824826-Z824826*l824826)- loads([38]) + diag(v826 * Cbus([38],[38])) == 0;

v828 == v824 - S824828*ctranspose(Z824828) - Z824828*ctranspose(S824828) + Z824828*l824828*ctranspose(Z824828);
[v824, S824828; ctranspose(S824828), l824828] >= 0;
diag(A *(S824828-Z824828*l824828) * AH)- loads([39, 40, 41]) + diag(A * v828 * Cbus([39, 40, 41],[39, 40, 41]) * AH) == diag(A * S828830 * AH) + 0;

v830 == v828 - S828830*ctranspose(Z828830) - Z828830*ctranspose(S828830) + Z828830*l828830*ctranspose(Z828830);
[v828, S828830; ctranspose(S828830), l828830] >= 0;
diag(A *(S828830-Z828830*l828830) * AH)- loads([42, 43, 44]) + diag(A * v830 * Cbus([42, 43, 44],[42, 43, 44]) * AH) == diag(A * S830854 * AH) + 0;

v854 == v830 - S830854*ctranspose(Z830854) - Z830854*ctranspose(S830854) + Z830854*l830854*ctranspose(Z830854);
[v830, S830854; ctranspose(S830854), l830854] >= 0;
diag(A *(S830854-Z830854*l830854) * AH)- loads([45, 46, 47]) + diag(A * v854 * Cbus([45, 46, 47],[45, 46, 47]) * AH) == [0; diag(S854856); 0] + diag(A * S854852 * AH) + 0;

v858 == v832 - S832858*ctranspose(Z832858) - Z832858*ctranspose(S832858) + Z832858*l832858*ctranspose(Z832858);
[v832, S832858; ctranspose(S832858), l832858] >= 0;
diag(A *(S832858-Z832858*l832858) * AH)- loads([51, 52, 53]) + diag(A * v858 * Cbus([51, 52, 53],[51, 52, 53]) * AH) == [diag(S858864); 0; 0] + diag(A * S858834 * AH) + 0;

v888 == Xfmr_ratio^2 * (v832 - S832888*ctranspose(Z832888) - Z832888*ctranspose(S832888) + Z832888*l832888*ctranspose(Z832888));
% v888 == v832 .* Xfmr_ratio^2;
[v832, S832888; ctranspose(S832888), l832888] >= 0;
diag(A *(S832888-Z832888*l832888) * AH)- loads([84, 85, 86]) + diag(A * v888 * Cbus([84, 85, 86],[84, 85, 86]) * AH) == diag(A * S888890 * AH) + 0;

v860 == v834 - S834860*ctranspose(Z834860) - Z834860*ctranspose(S834860) + Z834860*l834860*ctranspose(Z834860);
[v834, S834860; ctranspose(S834860), l834860] >= 0;
diag(A *(S834860-Z834860*l834860) * AH)- loads([57, 58, 59]) + diag(A * v860 * Cbus([57, 58, 59],[57, 58, 59]) * AH) == diag(A * S860836 * AH) + 0;

v842 == v834 - S834842*ctranspose(Z834842) - Z834842*ctranspose(S834842) + Z834842*l834842*ctranspose(Z834842);
[v834, S834842; ctranspose(S834842), l834842] >= 0;
diag(A *(S834842-Z834842*l834842) * AH)- loads([60, 61, 62]) + diag(A * v842 * Cbus([60, 61, 62],[60, 61, 62]) * AH) == diag(A * S842844 * AH) + 0;

v840 == v836 - S836840*ctranspose(Z836840) - Z836840*ctranspose(S836840) + Z836840*l836840*ctranspose(Z836840);
[v836, S836840; ctranspose(S836840), l836840] >= 0;
diag(A *(S836840-Z836840*l836840) * AH)- loads([66, 67, 68]) + diag(A * v840 * Cbus([66, 67, 68],[66, 67, 68]) * AH) == 0;

v862 == v836 - S836862*ctranspose(Z836862) - Z836862*ctranspose(S836862) + Z836862*l836862*ctranspose(Z836862);
[v836, S836862; ctranspose(S836862), l836862] >= 0;
diag(A *(S836862-Z836862*l836862) * AH)- loads([69, 70, 71]) + diag(A * v862 * Cbus([69, 70, 71],[69, 70, 71]) * AH) == [0; diag(S862838); 0] + 0;

v844 == v842 - S842844*ctranspose(Z842844) - Z842844*ctranspose(S842844) + Z842844*l842844*ctranspose(Z842844);
[v842, S842844; ctranspose(S842844), l842844] >= 0;
diag(A *(S842844-Z842844*l842844) * AH)- loads([72, 73, 74]) + diag(A * v844 * Cbus([72, 73, 74],[72, 73, 74]) * AH) == diag(A * S844846 * AH) + 0;

v846 == v844 - S844846*ctranspose(Z844846) - Z844846*ctranspose(S844846) + Z844846*l844846*ctranspose(Z844846);
[v844, S844846; ctranspose(S844846), l844846] >= 0;
diag(A *(S844846-Z844846*l844846) * AH)- loads([75, 76, 77]) + diag(A * v846 * Cbus([75, 76, 77],[75, 76, 77]) * AH) == diag(A * S846848 * AH) + 0;

v848 == v846 - S846848*ctranspose(Z846848) - Z846848*ctranspose(S846848) + Z846848*l846848*ctranspose(Z846848);
[v846, S846848; ctranspose(S846848), l846848] >= 0;
diag(A *(S846848-Z846848*l846848) * AH)- loads([78, 79, 80]) + diag(A * v848 * Cbus([78, 79, 80],[78, 79, 80]) * AH) == 0;

v816 == v850 - S850816*ctranspose(Z850816) - Z850816*ctranspose(S850816) + Z850816*l850816*ctranspose(Z850816);
[v850, S850816; ctranspose(S850816), l850816] >= 0;
diag(A *(S850816-Z850816*l850816) * AH)- loads([29, 31, 32]) + diag(A * v816 * Cbus([29, 31, 32],[29, 31, 32]) * AH) == [diag(S816818); 0; 0] + diag(A * S816824 * AH) + 0;

v832 == v852R - S852R832*ctranspose(Z852R832) - Z852R832*ctranspose(S852R832) + Z852R832*l852R832*ctranspose(Z852R832);
[v852R, S852R832; ctranspose(S852R832), l852R832] >= 0;
diag(A *(S852R832-Z852R832*l852R832) * AH)- loads([48, 49, 50]) + diag(A * v832 * Cbus([48, 49, 50],[48, 49, 50]) * AH) == diag(A * S832858 * AH) + diag(A * S832888 * AH) + 0;

v856 == v854_abc([2],[2]) - S854856*ctranspose(Z854856) - Z854856*ctranspose(S854856) + Z854856*l854856*ctranspose(Z854856);
[v854_abc([2],[2]), S854856; ctranspose(S854856), l854856] >= 0;
diag(S854856-Z854856*l854856)- loads([87]) + diag(v856 * Cbus([87],[87])) == 0;

v852 == v854 - S854852*ctranspose(Z854852) - Z854852*ctranspose(S854852) + Z854852*l854852*ctranspose(Z854852);
[v854, S854852; ctranspose(S854852), l854852] >= 0;
diag(A *(S854852-Z854852*l854852) * AH)- loads([88, 89, 90]) + diag(A * v852 * Cbus([88, 89, 90],[88, 89, 90]) * AH) == diag(A * S852852R * AH) + 0;

v864 == v858_abc([1],[1]) - S858864*ctranspose(Z858864) - Z858864*ctranspose(S858864) + Z858864*l858864*ctranspose(Z858864);
[v858_abc([1],[1]), S858864; ctranspose(S858864), l858864] >= 0;
diag(S858864-Z858864*l858864)- loads([91]) + diag(v864 * Cbus([91],[91])) == 0;

v834 == v858 - S858834*ctranspose(Z858834) - Z858834*ctranspose(S858834) + Z858834*l858834*ctranspose(Z858834);
[v858, S858834; ctranspose(S858834), l858834] >= 0;
diag(A *(S858834-Z858834*l858834) * AH)- loads([54, 55, 56]) + diag(A * v834 * Cbus([54, 55, 56],[54, 55, 56]) * AH) == diag(A * S834860 * AH) + diag(A * S834842 * AH) + 0;

v836 == v860 - S860836*ctranspose(Z860836) - Z860836*ctranspose(S860836) + Z860836*l860836*ctranspose(Z860836);
[v860, S860836; ctranspose(S860836), l860836] >= 0;
diag(A *(S860836-Z860836*l860836) * AH)- loads([63, 64, 65]) + diag(A * v836 * Cbus([63, 64, 65],[63, 64, 65]) * AH) == diag(A * S836840 * AH) + diag(A * S836862 * AH) + 0;

v838 == v862_abc([2],[2]) - S862838*ctranspose(Z862838) - Z862838*ctranspose(S862838) + Z862838*l862838*ctranspose(Z862838);
[v862_abc([2],[2]), S862838; ctranspose(S862838), l862838] >= 0;
diag(S862838-Z862838*l862838)- loads([92]) + diag(v838 * Cbus([92],[92])) == 0;

v890 == v888 - S888890*ctranspose(Z888890) - Z888890*ctranspose(S888890) + Z888890*l888890*ctranspose(Z888890);
[v888, S888890; ctranspose(S888890), l888890] >= 0;
diag(A *(S888890-Z888890*l888890) * AH)- loads([93, 94, 95]) + diag(A * v890 * Cbus([93, 94, 95],[93, 94, 95]) * AH) == 0;

v800 == vSOURCEBUS - SSOURCEBUS800*ctranspose(ZSOURCEBUS800) - ZSOURCEBUS800*ctranspose(SSOURCEBUS800) + ZSOURCEBUS800*lSOURCEBUS800*ctranspose(ZSOURCEBUS800);
[vSOURCEBUS, SSOURCEBUS800; ctranspose(SSOURCEBUS800), lSOURCEBUS800] >= 0;
diag(A *(SSOURCEBUS800-ZSOURCEBUS800*lSOURCEBUS800) * AH)- loads([4, 5, 6]) + diag(A * v800 * Cbus([4, 5, 6],[4, 5, 6]) * AH) == diag(A * S800802 * AH) + 0;

A * v814R * AH == (A * (v814([1, 2, 3],[1, 2, 3]) - S814814R*Z814814R' - Z814814R*S814814R' + Z814814R*l814814R*Z814814R') * AH) .* alphaM814R;
[v814([1, 2, 3],[1, 2, 3]), S814814R; ctranspose(S814814R), l814814R] >= 0;
diag(A *(S814814R-Z814814R*l814814R) * AH)- loads([23, 24, 25]) + diag(A * v814R * Cbus([23, 24, 25],[23, 24, 25]) * AH) == diag(A * S814R850 * AH) + 0;

A * v852R * AH == (A * (v852([1, 2, 3],[1, 2, 3]) - S852852R*Z852852R' - Z852852R*S852852R' + Z852852R*l852852R*Z852852R') * AH) .* alphaM852R;
[v852([1, 2, 3],[1, 2, 3]), S852852R; ctranspose(S852852R), l852852R] >= 0;
diag(A *(S852852R-Z852852R*l852852R) * AH)- loads([81, 82, 83]) + diag(A * v852R * Cbus([81, 82, 83],[81, 82, 83]) * AH) == diag(A * S852R832 * AH) + 0;

v808_abc == A * v808 * AH;
v862_abc == A * v862 * AH;
v824_abc == A * v824 * AH;
v854_abc == A * v854 * AH;
v816_abc == A * v816 * AH;
v858_abc == A * v858 * AH;


cvx_end


VSOURCEBUS = A * v0;
ISOURCEBUS800 = 1/trace(A * vSOURCEBUS([1, 2, 3],[1, 2, 3]) * AH) * ctranspose(A * SSOURCEBUS800 * AH)*VSOURCEBUS([1, 2, 3]);
V800 = VSOURCEBUS([1, 2, 3]) - A * ZSOURCEBUS800* AH *ISOURCEBUS800;
I800802 = 1/trace(A * v800([1, 2, 3],[1, 2, 3]) * AH) * ctranspose(A * S800802 * AH)*V800([1, 2, 3]);
V802 = V800([1, 2, 3]) - A * Z800802* AH *I800802;
I802806 = 1/trace(A * v802([1, 2, 3],[1, 2, 3]) * AH) * ctranspose(A * S802806 * AH)*V802([1, 2, 3]);
V806 = V802([1, 2, 3]) - A * Z802806* AH *I802806;
I806808 = 1/trace(A * v806([1, 2, 3],[1, 2, 3]) * AH) * ctranspose(A * S806808 * AH)*V806([1, 2, 3]);
V808 = V806([1, 2, 3]) - A * Z806808* AH *I806808;
I808810 = 1/trace(v808_abc([2],[2]) ) * ctranspose(S808810)*V808([2]);
V810 = V808([2]) - Z808810*I808810;
I808812 = 1/trace(A * v808([1, 2, 3],[1, 2, 3]) * AH) * ctranspose(A * S808812 * AH)*V808([1, 2, 3]);
V812 = V808([1, 2, 3]) - A * Z808812* AH *I808812;
I812814 = 1/trace(A * v812([1, 2, 3],[1, 2, 3]) * AH) * ctranspose(A * S812814 * AH)*V812([1, 2, 3]);
V814 = V812([1, 2, 3]) - A * Z812814* AH *I812814;
I814814R = 1/trace(A * v814([1, 2, 3],[1, 2, 3]) * AH) * ctranspose(A * S814814R * AH)*V814([1, 2, 3]);
V814R = (V814([1, 2, 3])  - A * Z814814R* AH *I814814R) .* alpha814R;
I814R850 = 1/trace(A * v814R([1, 2, 3],[1, 2, 3]) * AH) * ctranspose(A * S814R850 * AH)*V814R([1, 2, 3]);
V850 = V814R([1, 2, 3]) - A * Z814R850* AH *I814R850;
I850816 = 1/trace(A * v850([1, 2, 3],[1, 2, 3]) * AH) * ctranspose(A * S850816 * AH)*V850([1, 2, 3]);
V816 = V850([1, 2, 3]) - A * Z850816* AH *I850816;
I816818 = 1/trace(v816_abc([1],[1]) ) * ctranspose(S816818)*V816([1]);
V818 = V816([1]) - Z816818*I816818;
I816824 = 1/trace(A * v816([1, 2, 3],[1, 2, 3]) * AH) * ctranspose(A * S816824 * AH)*V816([1, 2, 3]);
V824 = V816([1, 2, 3]) - A * Z816824* AH *I816824;
I818820 = 1/trace(v818([1],[1]))*ctranspose(S818820)*V818([1]);
V820 = V818([1]) - Z818820*I818820;
I824826 = 1/trace(v824_abc([2],[2]) ) * ctranspose(S824826)*V824([2]);
V826 = V824([2]) - Z824826*I824826;
I824828 = 1/trace(A * v824([1, 2, 3],[1, 2, 3]) * AH) * ctranspose(A * S824828 * AH)*V824([1, 2, 3]);
V828 = V824([1, 2, 3]) - A * Z824828* AH *I824828;
I820822 = 1/trace(v820([1],[1]))*ctranspose(S820822)*V820([1]);
V822 = V820([1]) - Z820822*I820822;
I828830 = 1/trace(A * v828([1, 2, 3],[1, 2, 3]) * AH) * ctranspose(A * S828830 * AH)*V828([1, 2, 3]);
V830 = V828([1, 2, 3]) - A * Z828830* AH *I828830;
I830854 = 1/trace(A * v830([1, 2, 3],[1, 2, 3]) * AH) * ctranspose(A * S830854 * AH)*V830([1, 2, 3]);
V854 = V830([1, 2, 3]) - A * Z830854* AH *I830854;
I854856 = 1/trace(v854_abc([2],[2]) ) * ctranspose(S854856)*V854([2]);
V856 = V854([2]) - Z854856*I854856;
I854852 = 1/trace(A * v854([1, 2, 3],[1, 2, 3]) * AH) * ctranspose(A * S854852 * AH)*V854([1, 2, 3]);
V852 = V854([1, 2, 3]) - A * Z854852* AH *I854852;
I852852R = 1/trace(A * v852([1, 2, 3],[1, 2, 3]) * AH) * ctranspose(A * S852852R * AH)*V852([1, 2, 3]);
V852R = (V852([1, 2, 3]) - A * Z852852R* AH *I852852R) .* alpha852R;
I852R832 = 1/trace(A * v852R([1, 2, 3],[1, 2, 3]) * AH) * ctranspose(A * S852R832 * AH)*V852R([1, 2, 3]);
V832 = V852R([1, 2, 3]) - A * Z852R832* AH *I852R832;
I832858 = 1/trace(A * v832([1, 2, 3],[1, 2, 3]) * AH) * ctranspose(A * S832858 * AH)*V832([1, 2, 3]);
V858 = V832([1, 2, 3]) - A * Z832858* AH *I832858;
I832888 = 1/trace(A * v832([1, 2, 3],[1, 2, 3]) * AH) * ctranspose(A * S832888 * AH)*V832([1, 2, 3]);
V888 = (V832([1, 2, 3])  - A * Z832888* AH *I832888) * Xfmr_ratio;
% V888 = V832([1, 2, 3]) * Xfmr_ratio;
I858864 = 1/trace(v858_abc([1],[1]) ) * ctranspose(S858864)*V858([1]);
V864 = V858([1]) - Z858864*I858864;
I858834 = 1/trace(A * v858([1, 2, 3],[1, 2, 3]) * AH) * ctranspose(A * S858834 * AH)*V858([1, 2, 3]);
V834 = V858([1, 2, 3]) - A * Z858834* AH *I858834;
I888890 = 1/trace(A * v888([1, 2, 3],[1, 2, 3]) * AH) * ctranspose(A * S888890 * AH)*V888([1, 2, 3]);
V890 = V888([1, 2, 3]) - A * Z888890* AH *I888890;
I834860 = 1/trace(A * v834([1, 2, 3],[1, 2, 3]) * AH) * ctranspose(A * S834860 * AH)*V834([1, 2, 3]);
V860 = V834([1, 2, 3]) - A * Z834860* AH *I834860;
I834842 = 1/trace(A * v834([1, 2, 3],[1, 2, 3]) * AH) * ctranspose(A * S834842 * AH)*V834([1, 2, 3]);
V842 = V834([1, 2, 3]) - A * Z834842* AH *I834842;
I860836 = 1/trace(A * v860([1, 2, 3],[1, 2, 3]) * AH) * ctranspose(A * S860836 * AH)*V860([1, 2, 3]);
V836 = V860([1, 2, 3]) - A * Z860836* AH *I860836;
I842844 = 1/trace(A * v842([1, 2, 3],[1, 2, 3]) * AH) * ctranspose(A * S842844 * AH)*V842([1, 2, 3]);
V844 = V842([1, 2, 3]) - A * Z842844* AH *I842844;
I836840 = 1/trace(A * v836([1, 2, 3],[1, 2, 3]) * AH) * ctranspose(A * S836840 * AH)*V836([1, 2, 3]);
V840 = V836([1, 2, 3]) - A * Z836840* AH *I836840;
I836862 = 1/trace(A * v836([1, 2, 3],[1, 2, 3]) * AH) * ctranspose(A * S836862 * AH)*V836([1, 2, 3]);
V862 = V836([1, 2, 3]) - A * Z836862* AH *I836862;
I844846 = 1/trace(A * v844([1, 2, 3],[1, 2, 3]) * AH) * ctranspose(A * S844846 * AH)*V844([1, 2, 3]);
V846 = V844([1, 2, 3]) - A * Z844846* AH *I844846;
I862838 = 1/trace(v862_abc([2],[2]) ) * ctranspose(S862838)*V862([2]);
V838 = V862([2]) - Z862838*I862838;
I846848 = 1/trace(A * v846([1, 2, 3],[1, 2, 3]) * AH) * ctranspose(A * S846848 * AH)*V846([1, 2, 3]);
V848 = V846([1, 2, 3]) - A * Z846848* AH *I846848;


phasors=[];
phasors=[phasors;recover_voltage(VSOURCEBUS, 123)];
phasors=[phasors;recover_voltage(V800, 123)];
phasors=[phasors;recover_voltage(V802, 123)];
phasors=[phasors;recover_voltage(V806, 123)];
phasors=[phasors;recover_voltage(V808, 123)];
phasors=[phasors;recover_voltage(V810, 2)];
phasors=[phasors;recover_voltage(V812, 123)];
phasors=[phasors;recover_voltage(V814, 123)];
phasors=[phasors;recover_voltage(V814R, 123)];
phasors=[phasors;recover_voltage(V850, 123)];
phasors=[phasors;recover_voltage(V816, 123)];
phasors=[phasors;recover_voltage(V818, 1)];
phasors=[phasors;recover_voltage(V824, 123)];
phasors=[phasors;recover_voltage(V820, 1)];
phasors=[phasors;recover_voltage(V826, 2)];
phasors=[phasors;recover_voltage(V828, 123)];
phasors=[phasors;recover_voltage(V822, 1)];
phasors=[phasors;recover_voltage(V830, 123)];
phasors=[phasors;recover_voltage(V854, 123)];
phasors=[phasors;recover_voltage(V856, 2)];
phasors=[phasors;recover_voltage(V852, 123)];
phasors=[phasors;recover_voltage(V852R, 123)];
phasors=[phasors;recover_voltage(V832, 123)];
phasors=[phasors;recover_voltage(V858, 123)];
phasors=[phasors;recover_voltage(V888, 123) * Xfmr_ratio^-1];
phasors=[phasors;recover_voltage(V864, 1)];
phasors=[phasors;recover_voltage(V834, 123)];
phasors=[phasors;recover_voltage(V890, 123) * Xfmr_ratio^-1];
phasors=[phasors;recover_voltage(V860, 123)];
phasors=[phasors;recover_voltage(V842, 123)];
phasors=[phasors;recover_voltage(V836, 123)];
phasors=[phasors;recover_voltage(V844, 123)];
phasors=[phasors;recover_voltage(V840, 123)];
phasors=[phasors;recover_voltage(V862, 123)];
phasors=[phasors;recover_voltage(V846, 123)];
phasors=[phasors;recover_voltage(V838, 2)];
phasors=[phasors;recover_voltage(V848, 123)];

% change to per unit
phasors(:, 1) = phasors(:, 1) / Vbase;
phasors(:, 3) = phasors(:, 3) / Vbase;
phasors(:, 5) = phasors(:, 5) / Vbase;




Voltage_output=[];
Voltage_output = [Voltage_output; recover_voltage(V860, 123)];
Voltage_output = [Voltage_output; recover_voltage(V840, 123)];
Voltage_output = [Voltage_output; recover_voltage(V844, 123)];
Voltage_output = [Voltage_output; recover_voltage(V848, 123)];
Voltage_output = [Voltage_output; recover_voltage(V890* Xfmr_ratio^-1, 123)];
% Voltage_output = [Voltage_output; recover_voltage(V890, 123) ];
Voltage_output = [Voltage_output; recover_voltage(V830, 123)];
Voltage_output = [Voltage_output; recover_voltage(V802, 123)];
Voltage_output = [Voltage_output; recover_voltage(V808, 123)];
Voltage_output = [Voltage_output; recover_voltage(V818, 1)];
Voltage_output = [Voltage_output; recover_voltage(V820, 1)];
Voltage_output = [Voltage_output; recover_voltage(V816, 123)];
Voltage_output = [Voltage_output; recover_voltage(V824, 123)];
Voltage_output = [Voltage_output; recover_voltage(V824, 123)];
Voltage_output = [Voltage_output; recover_voltage(V828, 123)];
Voltage_output = [Voltage_output; recover_voltage(V854, 123)];
Voltage_output = [Voltage_output; recover_voltage(V832, 123)];
Voltage_output = [Voltage_output; recover_voltage(V858, 123)];
Voltage_output = [Voltage_output; recover_voltage(V858, 123)];
Voltage_output = [Voltage_output; recover_voltage(V834, 123)];
Voltage_output = [Voltage_output; recover_voltage(V860, 123)];
Voltage_output = [Voltage_output; recover_voltage(V836, 123)];
Voltage_output = [Voltage_output; recover_voltage(V862, 123)];
Voltage_output = [Voltage_output; recover_voltage(V842, 123)];
Voltage_output = [Voltage_output; recover_voltage(V844, 123)];
Voltage_output = [Voltage_output; recover_voltage(V846, 123)];
Voltage_output = [Voltage_output; recover_voltage(V806, 123)];
Voltage_output = [Voltage_output; recover_voltage(V810, 2)];
Voltage_output = [Voltage_output; recover_voltage(V820, 1)];
Voltage_output = [Voltage_output; recover_voltage(V822, 1)];
Voltage_output = [Voltage_output; recover_voltage(V824, 123)];
Voltage_output = [Voltage_output; recover_voltage(V826, 2)];
Voltage_output = [Voltage_output; recover_voltage(V828, 123)];
Voltage_output = [Voltage_output; recover_voltage(V830, 123)];
Voltage_output = [Voltage_output; recover_voltage(V856, 2)];
Voltage_output = [Voltage_output; recover_voltage(V858, 123)];
Voltage_output = [Voltage_output; recover_voltage(V864, 1)];
Voltage_output = [Voltage_output; recover_voltage(V834, 123)];
Voltage_output = [Voltage_output; recover_voltage(V860, 123)];
Voltage_output = [Voltage_output; recover_voltage(V836, 123)];
Voltage_output = [Voltage_output; recover_voltage(V840, 123)];
Voltage_output = [Voltage_output; recover_voltage(V838, 2)];
Voltage_output = [Voltage_output; recover_voltage(V844, 123)];
Voltage_output = [Voltage_output; recover_voltage(V846, 123)];
Voltage_output = [Voltage_output; recover_voltage(V848, 123)];

% change to per unit
Voltage_output(:, 1) = Voltage_output(:, 1) / Vbase;
Voltage_output(:, 3) = Voltage_output(:, 3) / Vbase;
Voltage_output(:, 5) = Voltage_output(:, 5) / Vbase;

% phasor890=recover_voltage(V890, 123) * 24.9/4.16;
% disp(V890);
% disp(phasor890)
disp(diag(A * S800802 * AH) / 1000);
% disp(diag(A * S824828 * AH) / 1000);
%disp(diag(A * S832858 * AH) / 1000);
%disp(diag(A * S832888 * AH) / 1000);
% disp(diag(S858864) / 1000);
% disp(diag(A * S842844 * AH) / 1000);
% disp(diag(A * S844846 * AH) / 1000);
% disp(diag(A * S846848 * AH) / 1000);
% disp(diag(S816818) / 1000);
% disp(diag(S818820) / 1000);
% disp(diag(S820822) / 1000);
% disp(loads([37]) / 1000);
% disp(eig([v820, S820822; ctranspose(S820822), l820822]));
% disp(diag(A * S816824 * AH) / 1000);
% disp(diag(A * S854852 * AH) / 1000);

% disp(eig([v800, S800802; ctranspose(S800802), l800802]));
% disp(eig([v802, S802806; ctranspose(S802806), l802806]));
% disp(eig([v806, S806808; ctranspose(S806808), l806808]));
% disp(eig([v808_abc([2],[2]), S808810; ctranspose(S808810), l808810]));
% disp(eig([v808, S808812; ctranspose(S808812), l808812]));
% disp(eig([v812, S812814; ctranspose(S812814), l812814]));
% disp(eig([v814R, S814R850; ctranspose(S814R850), l814R850]));
% disp(eig([v816_abc([1],[1]), S816818; ctranspose(S816818), l816818]));
% disp(eig([v816, S816824; ctranspose(S816824), l816824]));
% disp(eig([v818, S818820; ctranspose(S818820), l818820]));
% disp(eig([v820, S820822; ctranspose(S820822), l820822])); 
% disp(eig([v824_abc([2],[2]), S824826; ctranspose(S824826), l824826]));
% disp(eig([v824, S824828; ctranspose(S824828), l824828]));
% disp(eig([v828, S828830; ctranspose(S828830), l828830]));
% disp(eig([v830, S830854; ctranspose(S830854), l830854]));
% disp(eig([v832, S832858; ctranspose(S832858), l832858]));
% disp(eig([v832, S832888; ctranspose(S832888), l832888]));
% disp(eig([v834, S834860; ctranspose(S834860), l834860]));
% disp(eig([v834, S834842; ctranspose(S834842), l834842]));
% disp(eig([v836, S836840; ctranspose(S836840), l836840]));
% disp(eig([v836, S836862; ctranspose(S836862), l836862]));
% disp(eig([v842, S842844; ctranspose(S842844), l842844]));
% disp(eig([v844, S844846; ctranspose(S844846), l844846]));
% disp(eig([v846, S846848; ctranspose(S846848), l846848]));
% disp(eig([v850, S850816; ctranspose(S850816), l850816]));
% disp(eig([v852R, S852R832; ctranspose(S852R832), l852R832]));
% disp(eig([v854_abc([2],[2]), S854856; ctranspose(S854856), l854856]));
% disp(eig([v854, S854852; ctranspose(S854852), l854852]));
% disp(eig([v858_abc([1],[1]), S858864; ctranspose(S858864), l858864]));
% disp(eig([v858, S858834; ctranspose(S858834), l858834]));
% disp(eig([v860, S860836; ctranspose(S860836), l860836]));
% disp(eig([v862_abc([2],[2]), S862838; ctranspose(S862838), l862838]));
% disp(eig([v888, S888890; ctranspose(S888890), l888890]));