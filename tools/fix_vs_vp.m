function Model_out = fix_vs_vp(Model_in, Model_start)
% fixes vp and vs to the values from the starting model

            Model_vfix = change_parametrisation('rhomulambda','rhovsvp',Model_in);
            Modstart_vfix = change_parametrisation('rhomulambda','rhovsvp',Model_start);
            Model_vfix.vs = Modstart_vfix.vs;
            Model_vfix.vp = Modstart_vfix.vp;
            Model_out = change_parametrisation('rhovsvp','rhomulambda',Model_vfix);


end
