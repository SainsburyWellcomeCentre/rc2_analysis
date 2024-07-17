function animal_id = get_animal_id(probe_name)
%%GET_ANIMAL_ID Get the animal id from the probe name
%
%   ANIMAL_ID = get_animal_id(PROBE_NAME)
%   returns the correct animal ID from the probe namme.

animal_id = split(probe_name, "_");
if strcmp(probe_name, 'CA_176_1_rec1_rec2_rec3') || strcmp(probe_name, 'CA_176_3_rec1_rec2_rec3')

    animal_id = join([animal_id(1), animal_id(2), animal_id(3)], '_');
else
    animal_id = animal_id(1);
end