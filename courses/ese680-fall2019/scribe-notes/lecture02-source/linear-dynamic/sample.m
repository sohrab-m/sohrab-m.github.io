load iddata9 z9
nx = 5;
sys_arx = arx(z9, nx);
sys_n4sid = n4sid(z9, nx);

compare(z9,sys_arx,10)
compare(z9,sys_n4sid,10)