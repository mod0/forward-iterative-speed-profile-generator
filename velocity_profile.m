# Velocity Profile computation assuming Tracks are always available
l = 0.1; # step size = 0.1km
length_of_train = 1; # length of train in KM
track_lengths_in_km = [2,3,4,5];
permissible_speed_over_tracks = [20,10,5,30];
acceleration_in_tracks = [3,4,6,4];
deceleration_in_tracks = [-5,-3,-6,-2];
initial_train_velocity = 10;
final_train_velocity = 0;

number_of_steps_in_train_length = length_of_train / l;
number_of_steps_per_track = track_lengths_in_km / l;
number_of_steps = sum(number_of_steps_per_track);
number_of_steps_cumsum = cumsum(number_of_steps_per_track);
number_of_tracks = columns(track_lengths_in_km);

permissible_speed_array = zeros(1, number_of_steps);
max_permissible_speed_array = zeros(1, number_of_steps);
acceleration_across_steps_array = zeros(1,number_of_steps);
deceleration_across_steps_array = zeros(1,number_of_steps);


start_index = 1;
for i = 1:number_of_tracks
  permissible_speed_array(1,start_index:number_of_steps_cumsum(i)) = permissible_speed_over_tracks(i);
  # replace the arrays below with a dictionary or some other datastructure
  acceleration_across_steps_array(1,start_index:number_of_steps_cumsum(i)) = acceleration_in_tracks(i);
  deceleration_across_steps_array(1,start_index:number_of_steps_cumsum(i)) = deceleration_in_tracks(i);
  start_index = 1 + sum(number_of_steps_per_track(1:i));
endfor


# calculate maximum permissible speed array
for i = 1:number_of_steps
  start_index = i - number_of_steps_in_train_length;
  if(start_index < 1)
    start_index = 1;
  endif
  end_index = i;
  max_permissible_speed_array(i) = min(permissible_speed_array(1,start_index:end_index));
endfor

speed_profile_array = max_permissible_speed_array;

# Variables used in equations of motion
u = 0;
v = 0;

# final velocity - backward propagation
if(final_train_velocity < speed_profile_array(number_of_steps))
	v = final_train_velocity;
	u = sqrt(v^2 - 2 * deceleration_across_steps_array(number_of_steps) * l);
	if(u < speed_profile_array(number_of_steps))
		speed_profile_array(number_of_steps) = u;
	endif
endif


# backward propagation
for i = number_of_steps:-1:2
	if(speed_profile_array(i) < speed_profile_array(i-1))
		v = speed_profile_array(i);
		u = sqrt(v^2 - 2 * deceleration_across_steps_array(i-1) * l);
		if(u < speed_profile_array(i-1))
			speed_profile_array(i-1) = u;
		endif
	endif
endfor

# initial velocity - forward propagation
if(initial_train_velocity < speed_profile_array(1))
	u = initial_train_velocity;
	v = sqrt(u^2 + (2 * acceleration_across_steps_array(1) * l));
	if(v < speed_profile_array(1))
		speed_profile_array(1) = v;
	endif
endif

# forward propagation
for i = 1:(number_of_steps - 1)
	if(speed_profile_array(i) < speed_profile_array(i+1))
		u = speed_profile_array(i);
		v = sqrt(u^2 + (2 * acceleration_across_steps_array(i+1) * l));
		if(v < speed_profile_array(i+1))
			speed_profile_array(i+1) = v;
		endif
	endif
endfor


travel_time_per_step = l./speed_profile_array;
travel_time_cumulative = cumsum(travel_time_per_step);

subplot(5,1,1);
plot(permissible_speed_array);
xlabel('Distance Units');
ylabel('Speed');
subplot(5,1,2);
plot(max_permissible_speed_array);
xlabel('Distance Units');
ylabel('Speed');
subplot(5,1,3);
plot(permissible_speed_array - max_permissible_speed_array);
xlabel('Distance Units');
ylabel('Speed');
subplot(5,1,4);
plot(speed_profile_array);
xlabel('Distance Units');
ylabel('Speed');
subplot(5,1,5);
plot(travel_time_cumulative);
xlabel('Distance Units');
ylabel('CumulativeTime');
axis([1, 140 ,0, 20],'manual');