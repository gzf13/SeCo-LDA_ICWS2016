function writestreamstring4( WS , WO , ncharsperline , maxlines , fid , c , z , proba , focustopics , thr , NC );
  
L = [];
N = length( WS );
NL = 0;
for i=1:N
    if (WS( i ) == 0)
        L      = [ L '. ' ];
        spanst = 'normal';
        word   = '. ';
    else
        word   = WO{ WS( i )};
        status = c( i );
        topic  = z( i );
        prob   = proba( i );
        L      = [ L ' ' word ];
        spanst = 'filler';
        
        if (topic > 0)
            spanst = 'normal';
        end
        
        if (topic == focustopics)
           whcol = round( status * NC );
           spanst = sprintf( 't%d' , whcol );
        end
    end
    
    if length( L ) > ncharsperline
        NL = NL + 1;
        
        if (NL == maxlines+1)
            fprintf( fid , '...' );
            break
        else
            L = word;
            fprintf( fid , '\r\n' );
        end
    end

    fprintf( fid , ' <span class="%s">%s</span>' , spanst , word );
end

fprintf( fid , '\r\n' ); 

