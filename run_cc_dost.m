day = {'77','78','79','80','81','82','83','84','85','86'};
for iday=1:length(day)
    C1_A01_filt = eval(['C1_A01_',day{iday},'_filt']);
    cc_A05_A15_dost
    cc_05_15_filt_tmp(iday,:) = cc_05_15_filt;
    cc_05_15_tmp(iday,:) = cc_05_15;
end
