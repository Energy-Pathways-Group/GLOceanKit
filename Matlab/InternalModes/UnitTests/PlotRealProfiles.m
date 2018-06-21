load('SampleLatmixProfiles.mat');


for i=1:2:length(rhoProfile)
    im = InternalModes(rhoProfile{i},zProfile{i},zProfile{i},latitude);
    %     im.ShowLowestModesAtWavenumber(0);
    im.ShowLowestModesAtFrequency(im.f0);
end