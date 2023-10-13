function printButtonCallbackDFS(button, title, body)
        
    % remove button for printing
    set(button, 'Visible', 'off');
    drawnow;
    
    % body text default: nothing
    if ~exist('body','var'); body = ''; end

    % print figure to elog
    
    printfile = 'tempfile.png';
    print(printfile, '-dpng');
    printFigureToElogDFS(title, printfile, body);
    delete(printfile);

    % reset button
    set(button, 'Visible', 'on');
    
end