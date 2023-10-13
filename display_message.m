function display_message(destination,message,colour)

% display_message(message_location,message);

if nargin < 3
	colour = '';
end


if ishandle(destination)
    style = get(destination,'Style');
    switch style
        case 'text'
			if isempty(colour)
				set(destination,'String',message);
			else
				set(destination,'String',message,'BackgroundColor',colour);
			end
        case 'listbox'    
            old_message = get(destination,'String');
            new_message = old_message;
            new_message{end+1}=message;

            new_line_num = length(old_message)+1;
            set(destination,'String',new_message);
            set(destination,'Value',new_line_num);
    end
else
    disp(message);
end


end