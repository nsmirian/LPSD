
function printinlogbook(button, filename, Title, text )

text=['data seved in :', filename]
A=getframe(gcf);
imwrite(A.cdata, 'temporly.png');
result_log=hlc_send_to_logbook('title', Title, ...
    'text', text, 'image', 'temporly.png')

% delete('temporly.png')

end