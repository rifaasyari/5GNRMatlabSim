function SystemMain()
%SystemMain Run the expierment.
%To Do:
% *Learning on other spurs
% *Learning on main carrier.
% *Case for if CCs are different BWs
% *Add suppport for TX/RX on different freqs

%% Data Source
myLTE = LTE(5,'QPSK','uplink',200,-10e6);
myLTE.sampleArray = LTE.normalizeSignal(myLTE.sampleArray);

%% PA
myPA = PA(5);     %Set up a PA
%myDAC = DAC(10,9);
%myFrontend = Frontend(myPA,myDAC);
%myPA = WARP(1); %Set up WARP board

% Broadcast double pre signal
out = broadcast(myPA,myLTE.sampleArray);

%% DAC


%normal = LTE.normalizeSignal(out);
%LTE.plot_freqdomain(out,myLTE.CCs.CC1.systemFs,'','Post-PA');

%% DPD Processing
%Set up DPD unit and perform learning
%% Double
myDPD = SubBandDPD(myPA,myLTE,'IM3+');

%Apply learned DPD to signal
%signalWithDPD = applyDPDtoSignal(myDPD,myLTE);
DPDout1 = applyDPDtoSignal(myDPD,myLTE);

%Broadcast through PA
signalWithDPD = broadcast(myPA,DPDout1);

%normal = LTE.normalizeSignal(out);
LTE.plot_freqdomain(signalWithDPD,myLTE.CCs.CC1.systemFs,'','Double');


%% 8 bits
%  myDAC = DAC(8,7);
%  myFrontend = Frontend(myPA,myDAC);
%  myDPD = SubBandDPD(myFrontend,myLTE,'IM3+');
%  DPDout1 = applyDPDtoSignal(myDPD,myLTE);
%  signalWithDPD = broadcast(myPA,DPDout1);
%  LTE.plot_freqdomain(signalWithDPD,myLTE.CCs.CC1.systemFs,'','8 bits');
%  
%  %% 6 Bits
%  myDAC = DAC(6,5);
%  myFrontend = Frontend(myPA,myDAC);
%  myDPD = SubBandDPD(myFrontend,myLTE,'IM3+');
%  DPDout1 = applyDPDtoSignal(myDPD,myLTE);
%  signalWithDPD = broadcast(myPA,DPDout1);
%  LTE.plot_freqdomain(signalWithDPD,myLTE.CCs.CC1.systemFs,'','6 bits');
%  
%  %% 4 bits
%  myDAC = DAC(4,3);
%  myFrontend = Frontend(myPA,myDAC);
%  myDPD = SubBandDPD(myFrontend,myLTE,'IM3+');
%  DPDout1 = applyDPDtoSignal(myDPD,myLTE);
%  signalWithDPD = broadcast(myPA,DPDout1);
%  LTE.plot_freqdomain(signalWithDPD,myLTE.CCs.CC1.systemFs,'','4 bits');
%  
%  %% 2 bits
%  myDAC = DAC(2,1);
%  myFrontend = Frontend(myPA,myDAC);
%  myDPD = SubBandDPD(myFrontend,myLTE,'IM3+');
%  DPDout1 = applyDPDtoSignal(myDPD,myLTE);
%  signalWithDPD = broadcast(myPA,DPDout1);
%  LTE.plot_freqdomain(signalWithDPD,myLTE.CCs.CC1.systemFs,'','2 bits');

% 
%myDAC = DAC(1,0);
%myFrontend = Frontend(myPA,myDAC);
%myDPD = SubBandDPD(myPA,myLTE,'IM3+');
%DPDout1 = applyDPDtoSignal(myDPD,myLTE);
%signalWithDPD = broadcast(myPA,DPDout1);
%LTE.plot_freqdomain(signalWithDPD,myLTE.CCs.CC1.systemFs,'','1 bits');
%LTE.plot_freqdomain(out,myLTE.CCs.CC1.systemFs,'','No DPD');

%legend('show')

end

