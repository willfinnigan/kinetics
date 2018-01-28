import esterase_only.model_script
import Kinetics

model = esterase_only.model_script.model

y = model.run_model()

outputs = ["Ester", "Acid"]
Kinetics.print_model_output(y, model.time, model.species_names, outputs)
Kinetics.save_model_ouput(y, model, outputs, filename='')