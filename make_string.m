

str_txt = '.txt';
str_dash = ' - ';
str_slash = '\';
str0_d1 = 'd_1';
str0_d2 = 'd_2';
str0_d3 = 'd_3';
str0_R = 'R';
str0_Ax = 'Axes';

str0_h_Seq = 'hist_Seq';

folder_name = 'Data_1';

str1_beta = 'beta = ';
str2_beta=num2str(beta);
str3_beta=num2str(beta_seq);
str_beta = [str1_beta str2_beta str_dash  str3_beta];

str1_tail = 'tail = ';
str2_tail1 =num2str(tail1);
str2_tail2 =num2str(tail2);

str1_bmin = 'bond_min = ';
str2_bmin = num2str(bond_min);

str1_bmax = 'bond_max = ';
str2_bmax = num2str(bond_max);

str_force = 'force = ';
str_f1 = num2str(f1);
str_f2 = num2str(f2);
str_f3 = num2str(f3);

str1_Run = 'Run = ';
str2_Run = num2str(Run);

folder_path_0 = [folder_name str_slash str1_tail str2_tail1  str_dash str2_tail2 str_dash str_beta];
folder_path_2 = [str1_bmin str2_bmin str_dash str1_bmax str2_bmax];
folder_path_1 = [str_force str_f1 str_dash str_f2 str_dash str_f3];

folder_path = [folder_path_0 str_slash folder_path_1 str_slash folder_path_2];

str01 = [str1_bmin str2_bmin str_dash str1_bmax str2_bmax str_dash str_force str_f1 str_dash str_f2 str_dash str_f3 str_dash str1_tail str2_tail1 str_dash str2_tail2 str_dash str_beta str_dash str1_Run str2_Run];
str02 = [str1_bmin str2_bmin str_dash str1_bmax str2_bmax str_dash str_force str_f1 str_dash str_f2 str_dash str_f3 str_dash str1_tail str2_tail1 str_dash str2_tail2 str_dash str_beta];

str_d1 = [folder_path str_slash str0_d1 str_dash str01  str_txt];
str_d2 = [folder_path str_slash str0_d2 str_dash str01  str_txt];
str_d3 = [folder_path str_slash str0_d3 str_dash str01  str_txt];
str_R = [folder_path str_slash str0_R str_dash str01  str_txt];
str_Ax = [folder_path str_slash str0_Ax str_dash str01 str_txt];

str_h_Seq = [folder_path str_slash str0_h_Seq str_dash str01 str_txt];

