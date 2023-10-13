function [] = printFigureToElogDFS(title, image_filename, body)

    % define content
    author = '';
    keyword = 'Experiments';
    location = 'FLASH2';
    image = base64file(image_filename);
    
    % Create elog entry as XML  
    elogXML = createElogEntryXML(author, title, 'NONE', keyword, body, location, image);
    
    
    % '-l' is needed on Mac OS to enable raw print to a generic postscript queue RK 25.11.2013
    unix(['lpr -l -Pttflog ' elogXML]);
    
end

