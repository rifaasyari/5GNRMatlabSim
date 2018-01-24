classdef NR
   % NR class for generating physical layer signal
   % Based on:
   % 3GPP TS 38.201: "NR; Physical layer; General description"
   % 3GPP TS 38.202: "NR; Services provided by the physical layer"
   % 3GPP TS 38.211: "NR; Physical channels and modulation"
   %
   % RB = 12 subcarriers
   % RadioFrame = 10 ms / 10 subframes
   % Modualation from QPSK to 256 QAM
   % The channel coding scheme for transport blocks is quasi-cyclic LDPC codes with 2 base graphs and 8 sets of parity check matrices for each base graph, respectively.
   % The channel coding scheme for PBCH and control information is Polar coding based on nested sequences
   %
   % numerologyMu is the 3GPP parameter for subcarrier spacing etc
   % nRB can be between 24 and 275 for \mu = 0-3. for \mu = 4 the max is
   % 138 RBs and for \mu = 5 the max is 69. See:
   % < 38.211-Table 4.4.2-1: Minimum and maximum number of resource blocks.>
   %
   % MCS is the modulation and coding scheme. The modulation order is from
   % 38.214-Table5.1.3.1-2 index table 2 for PDSCH. We assume that the
   % higher layer parameter MCS-Table-PDSCH is set to '256QAM'
   
   properties (Constant)
      SUBCARRIERS_PER_RB = 12
      SYMBOLS_PER_SLOT   = 14
   end
   
   properties
      numerologyMu
      nRB
      BW % Bandwidth in MHz
      MCS
      modulationOrder
      subcarrierSpacing % spacing in kHz
      cyclicPrefixType
      slotsPerFrame
      slotsPerSubframe
      fftSize
      samplingRate % MHz
      nFrames
      nSymbols
      nSubCarrires
      PAPR
      signal
      cellId
   end
   
   methods
      function obj = NR(numerologyMu, nRB, MCS)
         if nargin == 0
            numerologyMu = 0;
            nRB = 100;
            MCS = 20;
         end
         
         obj = constructFrameStructure(obj, numerologyMu, nRB);
         obj = deconstructMCS(obj, MCS);
         
         obj.nFrames = 1;
         obj.cellId = 0;
         obj = signalGenerator(obj);
         
      end
      
      function obj = constructFrameStructure(obj,numerologyMu, nRB)
         obj.numerologyMu = numerologyMu;
         obj.nRB = nRB;
         obj.subcarrierSpacing = 2^numerologyMu * 15;
         obj.slotsPerFrame = 2^numerologyMu *10;
         obj.slotsPerSubframe = 2^numerologyMu;
         obj.nSymbols = obj.slotsPerFrame * obj.SYMBOLS_PER_SLOT
         if numerologyMu == 2
            obj.cyclicPrefixType = 'Extended';
         else
            obj.cyclicPrefixType = 'Normal';
         end
         obj.BW = nRB * obj.subcarrierSpacing * obj.SUBCARRIERS_PER_RB / 1000;
         obj.fftSize = getFFTSize(obj,nRB);
         obj.samplingRate = getSamplingRate(obj);
         obj.nSubCarrires = nRB * obj.SUBCARRIERS_PER_RB;
      end
      
      function obj = deconstructMCS(obj, MCS)
         load('38.214.Table5.1.3.1-2.mat');
         obj.MCS = MCS;
         obj.modulationOrder = MCSTable(MCS+1,2);
      end
      
      function bins = getFFTSize(obj,nRB)
         subcarriers = nRB * obj.SUBCARRIERS_PER_RB;
         bins = 2^(ceil(log2(subcarriers)));
      end
      
      function samplingRate = getSamplingRate(obj)
         samplingRate = obj.fftSize * obj.subcarrierSpacing / 1000;
      end
      
      function symbols = modulate(obj,preModulation)
         % From 38.211 Section 5.1
         symbols = ones(length(preModulation),1);
         if obj.modulationOrder == 2 % QPSK
            symbols = 2^(-0.5)*(...
               (1-2*preModulation(:,1)) + ...
               (1-2*preModulation(:,2))*1i);
         elseif obj.modulationOrder == 4 % 16QAM
            symbols = 10^(-0.5)*(...
               (1 - 2*preModulation(:,1)) .* (2 - (1-2*preModulation(:,3))) + ...
               (1 - 2*preModulation(:,2)) .* (2 - (1-2*preModulation(:,4)))*1i);
         elseif obj.modulationOrder == 6 % 64QAM
            symbols = 42^(-0.5)*(...
               (1-2*preModulation(:,1)).*...
               (4-(1-2*preModulation(:,3)).*(2-(1-2.*preModulation(:,5)))) + ...
               (1-2*preModulation(:,2)).*...
               (4-(1-2*preModulation(:,4)).*(2-(1-2.*preModulation(:,6))))*1i)  ;
            
         elseif obj.modulationOrder == 8 % 256QAM
            symbols = 170^(-0.5)*(...
               (1-2*preModulation(:,1)).*...
               (8-(1-2*preModulation(:,3)).*...
               (4-(1-2*preModulation(:,5)).*(2-(1-2.*preModulation(:,7))))) + ...
               (1-2*preModulation(:,2)).*...
               (8-(1-2*preModulation(:,4)).*...
               (4-(1-2*preModulation(:,6)).*(2-(1-2.*preModulation(:,8)))))*1i);
         else
            error('Invalid Modulation Order')
         end
      end
      
      function obj = signalGenerator(obj)
         nBits = 1*obj.modulationOrder;
         rawBits = randi([0 1],1,nBits);
         preModulation = reshape(rawBits,...
            [obj.modulationOrder,nBits/obj.modulationOrder])';
         symbolsToTransmit = modulate(obj,preModulation);
         %I'd prefer to change this to be actual resource elements and only
         %the number that there are. Then we'll flip things around for the
         %IFFT and zero pad. It'd make more sence and be easier to find
         %things.
         resourceElementGrid = zeros(obj.fftSize, obj.nSymbols);
         resourceElementGrid = generatePSS(obj,resourceElementGrid);
         
         %Fill the resource elements with the data in the queue
         
         %Fill first symbol.
         for j = 1:obj.nSymbols
            for i = 1:obj.nSubCarrires/2
               if resourceElementGrid(i+1,j) == 0
                  if isempty(symbolsToTransmit)
                     break;
                  end
                  resourceElementGrid(i+1,j) = symbolsToTransmit(1);
                  symbolsToTransmit = symbolsToTransmit(2:end);
                  if isempty(symbolsToTransmit)
                     break;
                  end
                  resourceElementGrid(end-i+1,j) = symbolsToTransmit(1);
                  symbolsToTransmit = symbolsToTransmit(2:end);
                  if isempty(symbolsToTransmit)
                     break;
                  end
               end
            end
         end
         
         resourceGrid = reshape(symbolsToTransmit,obj.nSubCarrires,[]);
         %resourceElementGrid(2:2+obj.nSubCarrires/2-1,:) = resourceGrid(1:obj.nSubCarrires/2,:);
         %resourceElementGrid(end - obj.nSubCarrires/2+1:end,:) = resourceGrid(obj.nSubCarrires/2+1:obj.nSubCarrires,:);
         
         
         %ifftInput(ceil((obj.fftSize - obj.nSubCarrires)/2)+1: ...
         %   floor((obj.fftSize - obj.nSubCarrires)/2) + obj.nSubCarrires+1,:)...
         %   = resourceGrid;
         ofdmData = ifft(resourceElementGrid)*sqrt(obj.nSubCarrires);
         ofdmData = insertCyclicPrefix(obj,ofdmData);
         obj.signal = ofdmData(:);
         
         spectrogram(obj.signal,kaiser(256,5),220,512,obj.samplingRate*10^6,'centered')
         %plotPSD(obj,obj.signal);
         %obj.PAPR = calculatePAPR(obj,obj.signal);
         %Maybe should put a CCDF plot
         %ccdf = comm.CCDF('PAPROutputPort',true,'MaximumPowerLimit', 50);
         %ccdf(obj.signal);
         %figure
         %plot(ccdf)
      end
      
      function out = insertCyclicPrefix(obj,signal)
         switch obj.cyclicPrefixType
            case 'Extended'
               nCyclicSamples = 512 * 2^-obj.numerologyMu;
               out = [signal; signal(end-nCyclicSamples-1:end,:)];
            case 'Normal'
               nCyclicSamples = 144 * 2^-obj.numerologyMu;
               out = [signal; signal(end-nCyclicSamples-1:end,:)];
         end
      end
      
      function out = generatePSS(obj,resourceGrid)
         % This is adapted from http://www.sharetechnote.com/html/Handbook_LTE_PSS.html
         u_shift = [25 29 34];
         NID = mod(obj.cellId,3);
         pssSequence = zeros(1,62);
         for n = 0:61
            u = u_shift(NID+1);
            if n <= 30
               d = exp(-j*pi*u*n*(n+1)/63);
            else
               d = exp(-j*pi*u*(n+1)*(n+2)/63);
            end
            pssSequence(n+1) = d;
         end
         
         % Insert the PSS into the resourceGrid.
         %
         % We assume that carrierfrequency <= 3 Ghz.
         % See 38.213 - 4.1 for more information
         % I'm not sure how often the PBCH block should be transmitted. I
         % just send on the first candidate symbol of a frame.
         if obj.subcarrierSpacing == 15  % 15 kHz subcarrier spacing
            pssSymbol = 2;
         elseif obj.subcarrierSpacing == 30
            pssSymbol = 4;
         elseif obj.subcarrierSpacing == 60
            pssSymbol = 2;
         else
            error('Cant do that subcarrier spacing')
         end
         
         % skip DC subcarrier. Put the last 32 of PSS here.
         resourceGrid(2:32,pssSymbol) = pssSequence(32:62);
         resourceGrid(end-30:end,pssSymbol) = pssSequence(1:31);
         out = resourceGrid;
      end
      
      function plotPSD(obj,signal)
         [pxx,f] = pwelch(signal,500,300,500,obj.samplingRate*10^6,'centered','power');
         figure
         plot(f/10^6,10*log10(pxx),'DisplayName','');
         xlabel('Frequency (MHz)')
         ylabel('Magnitude (dB)')
         grid on;
         hold on;
         axis([f(1)/10^6 f(end)/10^6 -inf inf])
      end
      
      function PAPR = calculatePAPR(~,signal)
         PAPR = 10*log10(max(abs(signal))^2/rms(signal)^2);
      end
   end
end