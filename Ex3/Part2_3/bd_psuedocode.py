MIN_POS = -0.5
MAX_POS = 4.5
FORCE_VECTOR = (-1,-0.5)
DIFFUSION_COEFFICENT = D


def BD_psuedocode(num_of_steps, delta_t, kbt):
	start_configuration = Vector(random(MIN_POS, MAX_POS),  random(MIN_POS, MAX_POS))
	current = start_configuration
	movement_list = []
	movement_list.append(current)
	
	for i in range(num_of_steps):
		current += delta_t / kbt * DIFFUSION_COEFFICENT * FORCE_VECTOR+ math.sqrt(6*DIFFUSION_COEFFICENT * delta_t) * random2DVector.Norm(0,1)
		current[0] = max(MIN_POS, current[0])
		current[0] = min(MAX_POS, current[0])
		
		current[1] = max(MIN_POS, current[1])
		current[1] = min(MAX_POS, current[1])
		movement_list.append(current)
	
	return movement_list