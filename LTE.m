classdef LTE
   % LTE Wrapper class to generate a combined LTE phy signal that'll add
   % multiple component carriers and shift them in baseband apprpriately.
   
   properties
      CCs
      direction
      nSymbols
      sampleArray
   end
   
   methods
      function obj = LTE(bw,modulation,direction,nSymbols,f)
         %LTE Sets up obj with global parameters for the entire signal and
         %builds the 1st component carrier.
         % direction: str, 'uplink' or 'downlink'
         % nSymbols: number of symbols to make
         % f: CenterFreq of 1st component carrier
         
         %These are global to everything in this LTE Signal
         obj.direction = direction;
         obj.nSymbols  = nSymbols;
         
         %Make first CC
         obj = makecomponentcarrier(obj,1,bw,modulation,f);
         n   = 1:length(obj.CCs.CC1.signalArray);
         obj.sampleArray = obj.CCs.CC1.signalArray...
            .*exp(2*pi*1j*n*f/obj.CCs.CC1.systemFs).';
      end
      
      function obj = newcomponentcarrier(obj,index,bw,modulation,f)
         %Create additional component carrier
         % f: CenterFreq of new component carrier
         obj = makecomponentcarrier(obj,index,bw,modulation,f);
         
         
         %ASSUMING SAME SAMPLE RATES!!!!
         
         %Need to check them. Compare and interpolate one with respect to
         %the other.
         
         indexStr = sprintf('CC%d',index);
         %Update the signal
         n   = 1:length(obj.CCs.(indexStr).signalArray);
         obj.sampleArray = obj.sampleArray + ...
            obj.CCs.(indexStr).signalArray...
            .*exp(2*pi*1j*n*f/obj.CCs.(indexStr).systemFs).';
      end
      
      function obj = makecomponentcarrier(obj,index,bw,modulationType,f)
         % purely makes a CC. Does not care about any other carriers or the
         % LTE signal as a whole.
         % inputs:
         %  obj   = the object.
         %  index = Index for this component carrier. 1 if 1st. 2 if
         %          2nd, etc.
         %  bw    = bandwidth of this carrier in MHz.
         % modulationTyp: str, 'QPSK'
         % f: Baseband center freq in Hz of new component carrier
         
         indexStr = sprintf('CC%d',index);
         
         obj.CCs.(indexStr).bandwidth = bw;
         obj.CCs.(indexStr).modulationType = modulationType;
         obj.CCs.(indexStr).centerFreq = f;
         
         [fftSize, nRB] = LTE.getfftsize(bw);
         SUBCARRIER_SPACING  = 15e3;  % 15 kHz bw per subcarrier
         SUBCARRIERS_PER_RB  = 12;
         obj.CCs.(indexStr).signalBandwidth = nRB * ...
            SUBCARRIER_SPACING*SUBCARRIERS_PER_RB;
         
         samplingFrequency = SUBCARRIER_SPACING * fftSize;
         nSubcarriers      = nRB * SUBCARRIERS_PER_RB;
         
         rrcRollOffFactor = 0.25; % Rolloff for pulse shaping
         
         % Upsampling for pulse shaping and possible Non-linearity modeling
         obj.CCs.(indexStr).upsamplingFactor = 21;%...
            %ceil(20*bw*1e6/samplingFrequency);
         
         %Define the sampling rate
         obj.CCs.(indexStr).systemFs = ...
            obj.CCs.(indexStr).upsamplingFactor*samplingFrequency;
         
         signal = LTE_SignalGenerator(obj,obj.CCs.(indexStr),nSubcarriers,...
            fftSize,rrcRollOffFactor);
         
         % Discard some samples from the end to make the signal cyclic
         obj.CCs.(indexStr).signalArray = ...
            signal(1:end - obj.CCs.(indexStr).upsamplingFactor*fftSize);
      end
      
      function LTE_signal = LTE_SignalGenerator(obj,CC,nSubcarriers, ...
            fftSize,rrcRollOffFactor)
         alphabet  = LTE.qamAlphabet(CC.modulationType); % Define symbols
         UpsamplingFactor = CC.upsamplingFactor;
         % Signal Generation
         switch obj.direction
            case 'downlink' % OFDM signal
               % symbol_RX                             =   alphabet(randi([1 length(alphabet)], users_bandwidth_subcarriers, obj.nSymbols));
               symbol_RX                               =   alphabet(ceil(length(alphabet)*rand(nSubcarriers,obj.nSymbols)));
               data_matrix_RX                          =   zeros(fftSize,obj.nSymbols);
               data_matrix_RX(2:2+nSubcarriers/2-1,:) = symbol_RX(1:1+nSubcarriers/2-1,:);
               data_matrix_RX(end-nSubcarriers/2+1:end,:) = symbol_RX(nSubcarriers/2+1:end,:);
               OFDM_data                               =   ifft(data_matrix_RX)*sqrt(nSubcarriers);
               signal_RX                               =   OFDM_data(:);
               pulse_RX                                =   rcosdesign(rrcRollOffFactor,50,UpsamplingFactor);
               s_over_RX                               =   zeros(UpsamplingFactor*length(signal_RX),1);
               s_over_RX(1:UpsamplingFactor:end)       =   signal_RX;
               S_IQ_RX                                 =   conv(pulse_RX,s_over_RX);
               LTE_signal                              =   S_IQ_RX((length(pulse_RX)+1)/2:end-(length(pulse_RX)+1)/2);
            case 'uplink' % SC-FDMA signal
               users_subchannel_locations              =   [-(nSubcarriers/2) (nSubcarriers/2)];
               users_subchannel_index                  =   transpose(users_subchannel_locations(1):users_subchannel_locations(2));
               users_N_subcarriers                     =   length(users_subchannel_index(:));
               % create QAM modulated symbols and put them in data_matrix
               % symbol_RX                             =   alphabet(randi([1 length(alphabet)], users_N_subcarriers, obj.nSymbols));
               symbol_RX                               =   alphabet(ceil(length(alphabet)*rand(users_N_subcarriers,obj.nSymbols)));
               symbol_dft                              =   fft(symbol_RX)/sqrt(length(symbol_RX));
               data_matrix                             =   symbol_dft;
               % Subcarrier Mapping
               index_vector                            =   mod(users_subchannel_index, fftSize) + 1;
               index_vector_blocks                     =   repmat(index_vector, 1, obj.nSymbols) + repmat(fftSize*(0:obj.nSymbols-1),  users_N_subcarriers, 1);
               symbol_mapped                           =   cell(fftSize,obj.nSymbols);
               load_zeros                              =   cellfun('isempty',symbol_mapped);
               symbol_mapped(load_zeros)               =   {0};
               symbol_mapped_double                    =   cell2mat(symbol_mapped);
               symbol_mapped_double(index_vector_blocks)=  data_matrix;
               % M-point IDFT
               s_ifft                                  =   sqrt(fftSize).*ifft(symbol_mapped_double);
               % P/S conversion
               signal                                  =   s_ifft(:);
               % Pulse shaping
               pulse                                   =   rcosdesign(rrcRollOffFactor,50,UpsamplingFactor); % rcosine(1,UpsamplingFactor,'sqrt',rrcRollOffFactor,50);
               s_over                                  =   zeros(UpsamplingFactor*length(signal),1);
               s_over(1:UpsamplingFactor:end)          =   signal;
               S_IQ                                    =   conv(pulse,s_over);
               LTE_signal                              =   S_IQ((length(pulse)+1)/2:end-(length(pulse)+1)/2);
         end
      end
   end
   methods(Static)
      function out = normalizeSignal(in)
         normalizedmax = 1;
         maxValue = max(max(real(in)),max(imag(in)));
         out = normalizedmax * maxValue^-1 .* in;
      end
      function alphabet = qamAlphabet(ModulationType)
         switch ModulationType
            case 'QPSK'
               MQAM = 4;
            case '16-QAM'
               MQAM = 16;
            case '64-QAM'
               MQAM = 64;
            otherwise
               error('Not Valid Modulation Scheme');
         end
         alphaMqam =   -(sqrt(MQAM)-1):2:(sqrt(MQAM)-1);
         A         =   repmat(alphaMqam,sqrt(MQAM),1);
         B         =   flipud(A');
         const_qam =   A+1j*B;
         const_qam =   const_qam(:);
         alphabet  =   const_qam;
      end
      
      function [fftSize, nRB] = getfftsize(bw)
         switch bw
            case 1.4
               fftSize = 128;
               nRB     = 6;
            case 3
               fftSize = 256;
               nRB     = 15;
            case 5
               fftSize = 512;
               nRB     = 25;
            case 10
               fftSize = 1024;
               nRB     = 50;
            case 15
               fftSize = 1536;
               nRB     = 75;
            case 20
               fftSize = 2048;
               nRB     = 100;
            otherwise
               error('Not Valid Bandwidth!');
         end
      end
      function f = plot_freqdomain(Rx_Signal,fs,colour,TITLE,POWER_PLOT_1MHZ,PA_Power_Measured)
         if nargin == 4
            POWER_PLOT_1MHZ = 0;
            PA_Power_Measured = 23;
         end
         if nargin == 1
            fs = 28800000;
            colour = 'k';
            TITLE = '';
            POWER_PLOT_1MHZ = 0;
            PA_Power_Measured = 23;
         end
         figure(100);
         
         if POWER_PLOT_1MHZ
            L = length(Rx_Signal);
            NFFT = L;%2^(nextpow2(L)-3); % Next power of 2 from length of Rx_Signal - 2 to reduce complexity
            Power_Rx_Signal = 10*log10(mean(abs(Rx_Signal).^2));
            PwrScaleFactor = 10^((PA_Power_Measured - Power_Rx_Signal)*0.1);
            Rx_Signal = Rx_Signal*sqrt(PwrScaleFactor);
            Rx_Signal_FD = fftshift(fft(Rx_Signal,NFFT)/L);
            Rx_Power_FD = (abs(Rx_Signal_FD).^2).';
            SubcarrierSpacing = fs/NFFT;
            % Bins_1MHz = round(UpsamplingFactor*1000/15);
            Bins_1MHz = round(1e6/SubcarrierSpacing);
            Rx_1MHz_Power = Bins_1MHz*smooth(Rx_Power_FD,Bins_1MHz,'moving');
            
            f = (fs/2)*linspace(-1,1,NFFT);
            % Plot double-sided amplitude spectrum.
            plot(f/1e6,10*log10(Rx_1MHz_Power),colour)
            title(TITLE)
            xlabel('Frequency (MHz)')
            ylabel('Power (dBm/MHz)')
            axis([-fs/(2*1e6) fs/(2*1e6) -Inf Inf])
            
         else
            %Nfft    = 2048;
            %Window  = kaiser(2000,9);
            %Signal_PSD = 10*log10(fftshift(pwelch(Rx_Signal,Window)));
            %plot((-1:2/Nfft:1-2/Nfft)*(fs/(2e6)),Signal_PSD,colour,'LineWidth',2);
            %xlabel('Frequency (MHz)')
            %ylabel('PSD')
            %axis([-fs/(2*1e6) fs/(2*1e6) -Inf Inf])
            Power_Rx_Signal = 10*log10(mean(abs(Rx_Signal).^2));
            PwrScaleFactor = 10^((PA_Power_Measured - Power_Rx_Signal)*0.1);
            Rx_Signal = Rx_Signal*sqrt(PwrScaleFactor);
            [pxx,f] = pwelch(Rx_Signal,500,300,500,fs,'centered','power');
            plot(f/10^6,10*log10(pxx),'DisplayName',TITLE);
            xlabel('Frequency (MHz)')
            ylabel('Magnitude (dB)')
            grid on;
            hold on;
            axis([f(1)/10^6 f(end)/10^6 -inf inf])
         end
      end
   end
end