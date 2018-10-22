function [simTime, tInt] = tensegSimTime(options,tEnd)

if (options.Refine)
    simTime = [0 tEnd];
    tInt = 0.01;
else
    prompt = 'Enter output time-step (e.g. 0.01): ';
    tInt = input(prompt);
    simTime = 0:tInt:tEnd;
end

end