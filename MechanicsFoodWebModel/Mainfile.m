%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Matlab Code supplementing the paper
% The mechanics of predator-prey interactions: first principles of physics predict predator-prey size ratios
% by Portalier, Fussmann, Loreau, Cherif
%
% Functional Ecology

% November 2018
%
% Matlab version: R2018b
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

%%%% READ ME %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The following code is used to run the model defined in the paper
% It calls several functions defined in other files
% It defines body masses for n species
% Then, the code computes for each species: 
% 1. A vector of traits (see fnSpecies function for details)
% The n vectors are stored into a cell array (SpeciesTraits)
% 2. An array of values for capture sequence (see fnMotion file for details)
% The n arrays are stored into a cell array (SpeciesArray)
%
% Several n*n matrices are computed to define all aspects of predator-prey relationship
% 1. Field metabolic rate (for predators)
% 2. Search cost
% 3. Capture cost
% 4. Handling cost
% 5. Energy content (for prey)
% 6. Net predation gains
% 7. Net predation gains per predator mass
%
% Details on calculation and full references for values are provided in 
% the main text or in the supplementary methods of the paper
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Physical parameters (sea water at 20 C) (change values for air, see main text)
BodyDensity=1080; % (kg.m-3)
MediumDensity=1026.95; % (kg.m-3)
DynamicViscosity=0.00139; % (N.s.m-2) 
Gravity=9.8; % (m.s-2)

param=[BodyDensity,MediumDensity,DynamicViscosity,Gravity];

%%%%%%%% Species body mass
% The example is done over 4 species only
% But any number of species should work
BodyMass=[1e-8,1e-5,1e-3,1e1];
NumberofSpecies=4;


%%%%%%%% Species traits
% see fnSpecies for details about the different values of the returned vector
SpeciesTraits=cell(1,NumberofSpecies);

for i=1:NumberofSpecies
    SpeciesTraits{i}=fnSpecies(BodyMass(i),param);
end

%%%%%%%% Array for capture sequence
SpeciesArray=cell(1,NumberofSpecies);

for i=1:NumberofSpecies
    % Retrieve species traits
    Sp=SpeciesTraits{i};
    % Time step for computation
    TimeStep=Sp(5)/100;
    % Switch variable for fnMotion function
    % i_Switch = 2: fnMotion is called for the computation of the capture sequence
    i_Switch=2;
    % Vector of parameters for fnMotion function (see fnMotion file for details)
    p=[Sp(1),BodyMass(i),Sp(3),Sp(4),Sp(5),Gravity,MediumDensity,DynamicViscosity,TimeStep,i_Switch];
    SpeciesArray{i}=fnMotion(BodyMass(i),p);
end
  
%%%%%%%% field metabolic rate (J.s-1)
fieldmetab=12.51*BodyMass'.^0.75';

%%%%%%%% Searching cost
% Defines a matrix for searching cost (output)
% and a intermediate matrix for encounter rate (used for searching time)  
searchingcost=zeros(NumberofSpecies,NumberofSpecies);
EncounterRate=zeros(NumberofSpecies,NumberofSpecies);

for i=1:NumberofSpecies
    Predator=SpeciesTraits{i};
    for j=1:NumberofSpecies
        Prey=SpeciesTraits{j};
        %%% prey abundance = 1% of space volume (see main text)
        Abundance=round(0.01/Prey(1));
        % If prey is very large, the minimal abundance should be 1 individual (not a fraction of individual)
        if Abundance<1.0
            Abundance=1;
        end
        %%% encounter rate
        % If the predator moves faster than the prey
        if Predator(10)>Prey(10)
            EncounterRate(i,j)=(Abundance*(pi*Predator(9)^2.0)*((Prey(10)^2.0)+3.0*(Predator(10)^2.0)))/(3.0*Predator(10));
        else
            % if the prey moves faster than the predator
            EncounterRate(i,j)=(Abundance*(pi*Predator(9)^2.0)*((Predator(10)^2.0)+3.0*(Prey(10)^2.0)))/(3.0*Prey(10));
        end
        %%% searching cost
        % cost = (motion cost at species-specific speed per time + metabolic rate per time) * searching time
        % searching time = 1 / encounter rate
        searchingcost(i,j)=(Predator(11)+fieldmetab(i))*(1/EncounterRate(i,j));
    end
end

%%%%%%%% Capture cost
% Defines a matrix for capture cost (output)
% and a intermediate matrix for capture probability 
% If prey is very large, the minimal abundance should be 1 individual (not
% capture cost = (capture cost per attempt + metabolic cost) * number of attempts
% number of attempts = 1 / capture probability
capturecost=zeros(NumberofSpecies,NumberofSpecies);
captureprobability=zeros(NumberofSpecies,NumberofSpecies);

for i=1:NumberofSpecies
    PredatorArray=SpeciesArray{i};
    Predator=SpeciesTraits{i};
    for j=1:NumberofSpecies
        PreyArray=SpeciesArray{j};
        Prey=SpeciesTraits{j};
        % time vectors
        predatortime=PredatorArray(4,:);
        preytime=PreyArray(4,:);
        % distance vectors
        predatordistance=PredatorArray(2,:);
        preydistance=PreyArray(2,:);
        % predator detection distance
        DetectiondistancePredator=Predator(9);
        % prey detection distance (initial distance when chase begins)
        DetectiondistancePrey=Prey(9);
        preydistance=preydistance+DetectiondistancePrey;
        % predator work
        predatorwork=PredatorArray(3,:);
            
        %%% find a crossing point
        % if Predator detection distance > Prey detection distance:
        % the predator can detect the prey before the prey can detect the predator 
        if DetectiondistancePredator>=DetectiondistancePrey
            % fit polynomial functions
            P1=polyfit(predatortime,predatordistance,5);
            P2=polyfit(preytime,preydistance,5);
            % check if coefficients are similar (i.e., parallel curves)
            i_check=abs(P1-P2);
            % parallel curves (no crossing point)
            if max(i_check(1:5))<1e-15
                distanceroot=NaN;
            else
                % if not parallel
                % determine the first root
                options=optimoptions(@fmincon,'Algorithm','active-set','MaxFunctionEvaluations',5000,'FunctionTolerance',1e-9);
                % constraints: min time >= 1e-10, max time <= predatortime(end) 
                A=[-1;1];
                b=[-1e-10;predatortime(end)];
                % starting point
                x0=1e-5;
                % Vector of parameters for fnVelocity function (see corresponding file for details)
                p=cell(1,2);
                p{1}=P1;
                p{2}=P2;
                % optimization
                root=fmincon(@(x)fnCapture(x,p),x0,A,b,[],[],[],[],[],options);
                % check whether or not curves cross (a root exists)
                if root<1e6
                     R1=polyval(P1,root);
                     R2=polyval(P2,root);
                     distanceroot=R2-R1;
                else
                    distanceroot=NaN;
                end
            end

            % cross or close enough for capture
            threshold=Predator(3)/100;
            if isnan(distanceroot)==0 && distanceroot < threshold
                % the predator can run over this period of time (feasible capture)
                if root<=predatortime(end)
                    % the prey can run over this period of time
                    if root<=preytime(end)
                        % find position in predator array
                        indextimepredator1=find(predatortime>=root);
                        indextimepredator=indextimepredator1(1);
                        % find position in prey array
                        indextimeprey1=find(preytime>=root);
                        indextimeprey=indextimeprey1(1);
                        capturecostperattempt=predatorwork(indextimepredator)+fieldmetab(i)*root;
                        % predator speed > 0
                        if PredatorArray(1,indextimepredator)>0
                            % capture probability
                            captureprobability(i,j)=1.0/(1.0+(PreyArray(1,indextimeprey)/PredatorArray(1,indextimepredator)));
                            % capture cost
                            capturecost(i,j)=capturecostperattempt/captureprobability(i,j);
                        else
                            % predator speed = 0
                            captureprobability(i,j)=0.0;
                            capturecost(i,j)=NaN;
                        end
                    else
                        % the prey stops before the predator can reach it
                        % prey speed = 0
                        captureprobability(i,j)=1.0;
                        capturecost(i,j)=capturecostperattempt;
                    end
                else
                    % the predator stops before
                    captureprobability(i,j)=0.0;
                    capturecost(i,j)=NaN;
                end
            else
                % no crossing point
                captureprobability(i,j)=0.0;
                capturecost(i,j)=NaN;
            end
        else
            % predator detection distance < prey detection distance
            % no crossing point: the prey escapes before detection
            captureprobability(i,j)=0.0;
            capturecost(i,j)=NaN;
        end
    end
end
   
%%%%%%%% Handling cost
% Defines a matrix for handling cost (output)
handlingcost=zeros(NumberofSpecies,NumberofSpecies);

for i=1:NumberofSpecies
    Predator=SpeciesTraits{i};
    for j=1:NumberofSpecies
        Prey=SpeciesTraits{j};
        %%%% handling time calculation
        %%% consumption time = number of bites * bite time
        % if predator bite size < prey size (multiple bites)
        if Predator(7)<Prey(2)
            consumptiontime=(Prey(2)/Predator(7))*Predator(8);
        else
            % prey consumed in one bite
            consumptiontime=Predator(8);
        end
        %%% digestion time
        digestiontime=2.3e4*(Prey(2)/Predator(2))*(Predator(2)^0.25);
        %%% handling time
        handlingtime=consumptiontime+digestiontime;
    
        %%%% handling cost calculation
        % Cumulated body volume
        Totalvolume=Predator(1)+Prey(1);
        % Cumulated body mass
        Totalmass=Predator(2)+Prey(2);
        % Max body radius
        Maxradius=max(Predator(3),Prey(3));
        % Max body section surface
        Maxbodysurface=max(Predator(4),Prey(4));
        % Time step for computation
        TimeStep=Predator(5)/100;
        % switch value for fnMotion function
        i_Switch=5;
        % Vector for fnHandlingMotion function (see corresponding file for details)
        p=[Totalvolume,Totalmass,Maxradius,Maxbodysurface,Predator(5),Gravity,MediumDensity,DynamicViscosity,TimeStep,i_Switch];
        % Work per time
        Workpt=fnMotion(Predator(2),p);
        %%% handling cost
        % function returns a cost (feasible interaction)
        if isnan(Workpt)~=1
            handlingcost(i,j)=(Workpt+fieldmetab(i))*handlingtime;
        else
            % function returns NA (non feasible interaction)
            handlingcost(i,j)=NaN;
        end
    end
end

%%%%%%%% Energy (J) given by the prey (in a matrix form)
energy=zeros(NumberofSpecies,NumberofSpecies);

for j=1:NumberofSpecies
    Prey=SpeciesTraits{j};
    energy(:,j)=Prey(6);
end

%%%%%%%% Calculation of net energetic gain for each predator-prey interaction
% Defines predation net gain
% and predation net gain per kg of predator
netgain=energy-(searchingcost+capturecost+handlingcost);

netgainperkg=zeros(NumberofSpecies,NumberofSpecies);
for i=1:NumberofSpecies
    Predator=SpeciesTraits{i};
    netgainperkg(i,:)=netgain(i,:)/Predator(2);
end
