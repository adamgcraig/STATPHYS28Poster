function [] = print_status_update_if_time( status_update_message, print_every_seconds_init )
%PRINT_STATUS_UPDATE_IF_TIME Print the message if enough time has elapsed.
%   On the first call, you can specify how often to print with
%   print_every_seconds_init.
%   Otherwise, it will print if 60 seconds have elapsed
%   since the last time it printed a message.

persistent timerVal code_start_time last_print_time print_every_seconds

if isempty(code_start_time)
    timerVal = tic;
    code_start_time = toc(timerVal);
    fprintf('setting up print_status_update_if_time at %g seconds\n', code_start_time)
end
if isempty(last_print_time)
    last_print_time = toc(timerVal);
end
if isempty(print_every_seconds)
    if ~exist('print_every_seconds_init', 'var')
        print_every_seconds_init = 60;
    end
    print_every_seconds = print_every_seconds_init;
end
current_time = toc(timerVal);
if current_time - last_print_time >= print_every_seconds
    disp( [sprintf('time %g, ', current_time) status_update_message] )
    last_print_time = current_time;
end

end