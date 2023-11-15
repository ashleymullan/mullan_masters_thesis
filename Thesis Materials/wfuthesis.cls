
%
% Wake Forest University - Thesis Document Class
% June 1998
% April 2008
% March 2011
% May 2012
%
% The WFU thesis document style file was modified to
% adhere to the new document class methods of LaTeX 2e.
% The following are a list of changes made:
%
%	macro definitions were changed to use \newcommand or
%	\renewcommand instead of \def whenever possible.
%	How to use \newcommand and \renewcommand:
%	\newcommand   {cmd}[args][opt]{def}
%	\renewcommand {cmd}[args][opt]{def}
%
%	  where 
%	  cmd  is the command name;
%	  args is an integer from 1 to 9 denoting the number of 
%	       arguments of the command being defined;
%	  opt  is included if the first argument is optional, and
%	       the first argument takes a default value of opt;
%	  def  is the text to be substituted for each occurance of cmd
%	       - #n indicates replacement with text of nth argument
%	       - can contain command or environment definitions
%	         only if the command being defined contains no arguments.
%
%	\def's were not replaced with \newcommand or \renewcommand
%	in the following cases:
%	  if the command name began with \end, or
%	  if the command definition contained command or environment
%	     definitions which contained arguments.
%	  if the command definition was preceded by \long
%
%	The commands \footheight, \@maxsep and \@dblmaxsep were
%	removed in LaTeX 2e, and they should no longer
%	be set in document classes.  These commands were removed from
%	the style document to create this class.
%
%	\rm, \sf, \tt, \bf, \it, \sl, and \sc are not defined in
%	the LaTeX 2e kernel, so they need to be defined by the
%	document class files.  They were defined in this class using
%	the New Font Selection Style so that they are orthogonal.
%	Also, they were defined so that the math version of these
%	commands will be automatically selected in math mode.	
%
%	Several fonts (vpt, vipt, etc.) are no longer preloaded in
%	LaTeX 2e, so these were defined in the document class. 		

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{wfuthesis}[2011/04/01 WFU Thesis class]

%--  draft option
%--
%-- \overfullrule = 0pt                % Default is don't mark overfull hboxes.
\newcommand{\textspace}{\openup 3\jot\relax}  % Default is double spaced text
\newcommand{\versiontag}{}                    % Default is NO tag


\DeclareOption{draft}{
	% -- DRAFT is double spaced
	% -- \newcommand\textspace{}  % Single spaced text
	\overfullrule 5pt 	% Causes overfull hboxes to be marked.
	\renewcommand\versiontag{{\normalsize\rm -
		 DRAFT \number\month/\number\day{} -}}
}

% Define three switches: \ifAM, \ifMS, and \ifPhD.  At most one must be
% true when the title page is actually printed.  They are all set false
% here, and one is normally set true by specifying an optional 
% \documentstyle parameter, one of: AM, MS, PhD.  If they are all false
% when \maketitle is invoked, the effect the same as \AMtrue.
\newif\ifAM  \AMfalse
\newif\ifMS  \MSfalse
\newif\ifPhD \PhDfalse

\DeclareOption{MS}{\MStrue}
\DeclareOption{AM}{\AMtrue}
\DeclareOption{MA}{\AMtrue}
\DeclareOption{PhD}{\PhDtrue}
\DeclareOption{Phd}{\PhDtrue}

%
% The \ProcessOptions command causes the execution of every FOO
% which is declared (\DeclarOption)  and for which the user typed the 
% FOO option in his
% \documentclass command.  For every option BAR he typed for which
% BAR is not defined, the file BAR.sty will be read after the present
% (main) .STY file is executed.
%

\ProcessOptions

%  ****************************************
%  *               FONTS                  *
%  ****************************************
%

\lineskip 1pt            	% \lineskip is 1pt for all font sizes.
\normallineskip 1pt
\renewcommand{\baselinestretch}{1}

%
% Each size-changing command \SIZE executes the command
%        \@setsize\SIZE{BASELINESKIP}\FONTSIZE\@FONTSIZE
% where:
%   BASELINESKIP = Normal value of \baselineskip for that size.  (Actual
%                  value will be \baselinestretch * BASELINESKIP.)
%
%  \FONTSIZE     = Name of font-size command.  The currently available
%                  (preloaded) font sizes are: \vpt (5pt), \vipt (6pt),
%                  \viipt (etc.), \viiipt, \ixpt, \xpt, \xipt, \xiipt,
%                  \xivpt, \xviipt, \xxpt, \xxvpt.
%		   THESE FONTS ARE NO LONGER PRELOADED! in LateX 2e
%
%  \@FONTSIZE    = The same as the font-size command except with an
%                  '@' in front---e.g., if \FONTSIZE = \xivpt then
%                  \@FONTSIZE = \@xivpt.
%

\newcommand{\vpt}{\fontsize{5}{6}\normalfont\selectfont}
\newcommand{\vipt}{\fontsize{6}{7}\normalfont\selectfont}
\newcommand{\viipt}{\fontsize{7}{8}\normalfont\selectfont}
\newcommand{\viiipt}{\fontsize{8}{9}\normalfont\selectfont}
\newcommand{\ixpt}{\fontsize{9}{10}\normalfont\selectfont}
\newcommand{\xpt}{\fontsize{10}{12}\normalfont\selectfont}
\newcommand{\xipt}{\fontsize{11}{13}\normalfont\selectfont}
\newcommand{\xiipt}{\fontsize{12}{14}\normalfont\selectfont}
\newcommand{\xivpt}{\fontsize{14}{16}\normalfont\selectfont}
\newcommand{\xviipt}{\fontsize{17}{20}\normalfont\selectfont}
\newcommand{\xxpt}{\fontsize{20}{24}\normalfont\selectfont}
\newcommand{\xxvpt}{\fontsize{25}{30}\normalfont\selectfont}

%
% Previously \@normalsize was defined instead of \normalsize, for
% reasons of efficiency. However, in LaTeX 2e \@normalsize cannot
% be used in this way, because \@normalsize is always redefined
% to the current definition of \normalsize.  
% \normalsize MUST be redefined in a class definition to be used.
%
\renewcommand{\normalsize}{\@setsize\normalsize{14.5pt}\xiipt\@xiipt
	\abovedisplayskip 12pt plus3pt minus7pt%
	\belowdisplayskip \abovedisplayskip
	\abovedisplayshortskip  \z@ plus3pt%
	\belowdisplayshortskip  6.5pt plus3.5pt minus3pt}

\newcommand{\small}{\@setsize\small{13.6pt}\xipt\@xipt
	\abovedisplayskip 11pt plus3pt minus6pt%
	\belowdisplayskip \abovedisplayskip
	\abovedisplayshortskip  \z@ plus3pt%
	\belowdisplayshortskip  6.5pt plus3.5pt minus3pt
	\renewcommand{\@listi}{\parsep 4.5pt plus 2pt minus 1pt
        \itemsep \parsep
        \topsep 9pt plus 3pt minus 5pt}}

\newcommand{\footnotesize}{\@setsize\footnotesize{12pt}\xpt\@xpt
	\abovedisplayskip 10pt plus2pt minus5pt%
	\belowdisplayskip \abovedisplayskip
	\abovedisplayshortskip  \z@ plus3pt%
	\belowdisplayshortskip  6pt plus3pt minus3pt
	\renewcommand\@listi{\topsep 6pt plus 2pt minus 2pt
	\parsep 3pt plus 2pt minus 1pt\itemsep \parsep}}

%
% \bf, \rm, etc. are no longer part of the LaTeX kernel so they need to
% be defined. These are defined here so they are orthogonal.
%
\DeclareOldFontCommand{\bf}{\bfseries}{\mathbf}
\DeclareOldFontCommand{\rm}{\rmfamily}{\mathrm}
\DeclareOldFontCommand{\sf}{\sffamily}{\mathsf}
\DeclareOldFontCommand{\tt}{\ttfamily}{\mathtt}
\DeclareOldFontCommand{\it}{\itshape}{\mathit}
\DeclareOldFontCommand{\sl}{\slshape}{\mathsl}
\DeclareOldFontCommand{\sc}{\scshape}{\mathsc}
%

\newcommand\scriptsize{\@setsize\scriptsize{9.5pt}\viiipt\@viiipt}
\newcommand\tiny{\@setsize\tiny{7pt}\vipt\@vipt}
\newcommand\large{\@setsize\large{18pt}\xivpt\@xivpt}
\newcommand\Large{\@setsize\Large{22pt}\xviipt\@xviipt}
\newcommand\LARGE{\@setsize\LARGE{25pt}\xxpt\@xxpt}
\newcommand\huge{\@setsize\huge{30pt}\xxvpt\@xxvpt}
\let\Huge=\huge


\normalsize  % Choose the normalsize font.

%
%  ****************************************
%  *            PAGE LAYOUT               *
%  ****************************************
%
% All margin dimensions measured from a point one inch from top and side
% of page.  Dimensions shrink by about 2 percent
%

%-- SIDE MARGINS:
\oddsidemargin 36pt    %  Left margin on odd-numbered pages.
\evensidemargin 36pt   %  Note that \oddsidemargin = \evensidemargin
\marginparwidth 68pt   %  Width of marginal notes.
\marginparsep 10pt     % Horizontal space between outer margin and marginal note


%-- VERTICAL SPACING: Modified 24 Nov 87 CRM
%-- \topmargin 0pt vice 20pt
                         % Top of page:
\topmargin 0pt           %    Nominal distance from top of page to top of
                         %    box containing running head.
%\headheight 12pt         %    Height of box containing running head.
\headheight 0pt
%\headsep 25pt            %    Space between running head and text.
\headsep 0pt
% \topskip = 10pt        %    '\baselineskip' for first line of page.
                         % Bottom of page:
\footskip 30pt           %    Distance from baseline of box containing foot
                         %    to baseline of last line of text.


%-- DIMENSION OF TEXT:
%--  8 Sept 87 -- changed the page sizes for new style req't
%--  crm@duke.
%----changed text height from 615 to 600 because text was 
%    sometimes running into page numbers.

\textheight 615pt      % Height of text (including footnotes and figures,
                         % excluding running head and foot). 635 vice
						 % 605
\textwidth 434pt         % Width of text line. 400 vice 430
                         % For two-column mode:
\columnsep 10pt          %    Space between columns
\columnseprule 0pt       %    Width of rule between columns.

%-- A \raggedbottom command causes 'ragged bottom' pages: pages set to
%-- natural height instead of being stretched to exactly \textheight.

%-- FOOTNOTES:

\footnotesep 14pt    % Height of strut placed at the beginning of every
                     % footnote.
\skip\footins 8pt plus 2pt minus 4pt  % Space between last line of text and
                                      % top of first footnote.

% 
% FLOATS: (a float is something like a figure or table)
% 
%  FOR FLOATS ON A TEXT PAGE:
%
%
\floatsep 14pt plus 2pt minus 4pt        % Space between adjacent floats moved
                                         % to top or bottom of text page.
\textfloatsep 20pt plus 6pt minus 4pt    % Space between main text and floats
                                         % at top or bottom of page.
\intextsep 14pt plus 4pt minus 4pt       % Space between in-text figures and
                                         % text.

%--  FOR FLOATS ON A SEPARATE FLOAT PAGE:
\@fptop 0pt plus 1fil    % Stretch at top of float page. (Must be 0pt plus ...)
\@fpsep 10pt plus 2fil   % Space between floats on float page.
\@fpbot 0pt plus 1fil    % Stretch at bottom of float page. (Must be 0pt plus ... )

%-- MARGINAL NOTES:
%
\marginparpush 7pt       % Minimum vertical separation between two marginal
                         % notes.

%
%  ****************************************
%  *           PARAGRAPHING               *
%  ****************************************
%
%
\parskip 2pt plus 1pt minus 1pt    % Extra vertical space between paragraphs.
\parindent 1.5em                   % Width of paragraph indentation.
\topsep 10pt plus 4pt minus 6pt    % Extra vertical space, in addition to
                                   % \parskip, added above and below list and
                                   % paragraphing environments.
\partopsep 3pt plus 2pt minus 2pt  % Extra vertical space, in addition to
                                   % \parskip and \topsep, added when user
                                   % leaves blank line before environment.
\itemsep 5pt plus 2.5pt minus 1pt  % Extra vertical space, in addition to
                                   % \parskip, added between list items.

%-- The following page-breaking penalties are defined

\@lowpenalty   51      % Produced by \nopagebreak[1] or \nolinebreak[1]
\@medpenalty  151      % Produced by \nopagebreak[2] or \nolinebreak[2]
\@highpenalty 301      % Produced by \nopagebreak[3] or \nolinebreak[3]

\@beginparpenalty -\@lowpenalty    % Before a list or paragraph environment.
\@endparpenalty   -\@lowpenalty    % After a list or paragraph environment.
\@itempenalty     -\@lowpenalty    % Between list items.


%    ****************************************
%    *        CHAPTERS AND SECTIONS         *
%    ****************************************
%
% DEFINE COUNTERS:
%
% \newcounter{NEWCTR}[OLDCTR] : Defines NEWCTR to be a counter, which is
%                               reset to zero when counter OLDCTR is stepped.
%                               Counter OLDCTR must already be defined.

\newcounter {chapter}
\newcounter {section}[chapter]
\newcounter {subsection}[section]
\newcounter {subsubsection}[subsection]
\newcounter {paragraph}[subsubsection]
\newcounter {subparagraph}[paragraph]

% For any counter CTR, \theCTR is a macro that defines the printed version
% of counter CTR.  It is defined in terms of the following macros:
%
%  \arabic{COUNTER} : The value of COUNTER printed as an arabic numeral.
%  \roman{COUNTER}  : Its value printed as a lower-case roman numberal.
%  \Roman{COUNTER}  : Its value printed as an upper-case roman numberal.
%  \alph{COUNTER}   : Value of COUNTER printed as a lower-case letter:
%                         1 = a, 2 = b, etc.
%  \Alph{COUNTER}   : Value of COUNTER printed as an upper-case letter:
%                           1 = A, 2 = B, etc.
%

\renewcommand{\thechapter}       {\arabic{chapter}}
\renewcommand{\thesection}       {\thechapter.\arabic{section}}
\renewcommand{\thesubsection}    {\thesection.\arabic{subsection}}
\renewcommand{\thesubsubsection} {\thesubsection .\arabic{subsubsection}}
\renewcommand{\theparagraph}     {\thesubsubsection.\arabic{paragraph}}
\renewcommand{\thesubparagraph}  {\theparagraph.\arabic{subparagraph}}


% \@chapapp is initially defined to be 'Chapter'.  The \appendix
% command redefines it to be 'Appendix'.

\newcommand{\@chapapp}{Chapter}

% \@makechapterhead {TEXT} : Makes the heading for the \chapter command.
%

\newcommand{\@makechapterhead}[1]{     % Heading for \chapter command
  	{ \parindent 0pt \raggedright
    	    \ifnum \c@secnumdepth >\m@ne     % IF secnumdepth > -1 THEN
      	    \large\bf \@chapapp{}\ \thechapter :  % Print 'Chapter' and number.
      	      \hskip 20pt	       %space between number and title
	    \fi                  
    	    \large \bf                 % Title.
    	    #1\par
    	    \nobreak                   % TeX penalty to prevent page break.
    	    \vskip 40pt                % Space between title and text.
  	}
  	%\advance\@colht-\footskip      % tell LaTeX to build a shorter page!
  	\pagegoal\@colht               % let TeX in on what's happening
}

% \@makeschapterhead {TEXT} : Makes the heading for the \chapter* command.
%
% 3-26-91 djj
% changed to center the heading at the top of the page

\newcommand{\@makeschapterhead}[1]{    % Heading for \chapter* command
	\begin{center}
	\large \bf
	#1
	\end{center}
	\vskip 40pt
}


% \secdef{UNSTARCMDS}{STARCMDS} :
%    When defining a \chapter or \section command without using
%    \@startsection, you can use \secdef as follows:
%       \def\chapter { ... \secdef \CMDA \CMDB }
%       \def\CMDA    [1][]{ ...   % Command to define \chapter[...]{...}
%       \def\CMDB    [1]{ ...       % Command to define \chapter*{...}

\newcommand{\chapter}{\clearpage      	% Starts new page.
   	\thispagestyle{headings}  % Page style of chapter page is 'headings'
   	\global\@topnum\z@        % Prevents figures from going at top of page.
   	\@afterindenttrue         % Indent the first paragraph.  Change
   	\secdef\@chapter\@schapter}   	% to \@afterindentfalse to
					% suppress indent.


\newcommand\@chapter[2][]{
	\ifnum \c@secnumdepth >\m@ne
            \refstepcounter{chapter}
            \typeout{\@chapapp\ \space\the\c@chapter.}
            \addcontentsline{toc}{chapter}{\@chapapp\ \protect
            \numberline{\thechapter} #1}
	\else
      	    \addcontentsline{toc}{chapter}{\@chapapp\ #1:}
	\fi
   	\chaptermark{#1}
   	\addtocontents{lof}{\protect\addvspace{10pt}} 
					% Adds between-chapter space
   	\addtocontents{lot}{\protect\addvspace{10pt}} 
					% to lists of figs & tables.
   	\@makechapterhead{#2}
   	\@afterheading             % Routine called after chapter and
}                                  % section heading.

\newcommand{\@schapter}[1]{\@makeschapterhead{#1} \@afterheading}

\newcommand{\pchapter}[1]{%% A preliminary chapter
	\chapter*{#1}
	\addcontentsline{toc}{chapter}{#1}
	\thispagestyle{headings}
	\typeout{#1}
}

% \@startsection {NAME}{LEVEL}{INDENT}{BEFORESKIP}{AFTERSKIP}{STYLE}
%            optional * [ALTHEADING]{HEADING}
%    Generic command to start a section.
%    NAME       : e.g., 'subsection'
%    LEVEL      : a number, denoting depth of section -- e.g., chapter=1,
%                 section = 2, etc.  A section number will be printed if
%                 and only if LEVEL < or = the value of the secnumdepth
%                 counter.
%    INDENT     : Indentation of heading from left margin
%    BEFORESKIP : Absolute value = skip to leave above the heading.
%                 If negative, then paragraph indent of text following
%                 heading is suppressed.
%    AFTERSKIP  : if positive, then skip to leave below heading,
%                       else - skip to leave to right of run-in heading.
%    STYLE      : commands to set style
%  If '*' missing, then increments the counter.  If it is present, then
%  there should be no [ALTHEADING] argument.  A sectioning command
%  is normally defined to \@startsection + its first six arguments.

% NOTE:  \section changed to \large (May, 2012)
%        \subsection changed to \normalsize (May, 2012)

\newcommand{\section}{\@startsection {section}{1}{\z@}
	{-3.5ex plus -1ex minus -.2ex}{2.3ex plus .2ex}{\large\bf}}
\newcommand{\subsection}{\@startsection{subsection}{2}{\z@}
	{-3.25ex plus -1ex minus -.2ex}{1.5ex plus .2ex}{\normalsize\bf}}
\newcommand{\subsubsection}{\@startsection{subsubsection}{3}{\z@}
	{-3.25ex plus -1ex minus -.2ex}{1.5ex plus .2ex}{\normalsize\bf}}
\newcommand{\paragraph}{\@startsection
     	{paragraph}{4}{\z@}{3.25ex plus 1ex minus .2ex}{-1em}{\normalsize\bf}}
\newcommand{\subparagraph}{\@startsection
     	{subparagraph}{4}{\parindent}{3.25ex plus 1ex minus
     	.2ex}{-1em}{\normalsize\bf}}

% Default initializations of \...mark commands.  (See below for their
% uses in defining page styles.
%

\newcommand{\chaptermark}[1]{}

% The value of the counter secnumdepth gives the depth of the
% highest-level sectioning command that is to produce section numbers.
%

\setcounter{secnumdepth}{2}

% COPYRIGHT
%
%
\newcommand{\Copyright}{
	\typeout{Copyright}
	\newpage\thispagestyle{empty}\null\vfil
	\vskip 0.75\textheight
	\begin{center}
		\def\thanks##1{} % temporarily kill any thanks lines
		{\large Copyright \copyright\ \number\year\ by \@author\par}
		{\large All rights reserved \par}
	\end{center}
	\par\vfil\null\newpage
}

% ABSTRACT
%
%   - this produces pages suitable for a UW-Abstract
%
\newcommand{\abstract}{\pchapter{Abstract}}

%  UNIV. Microfilms abstract
%
%
\newcommand{\UMIabstract}{
	\titlepage\pagestyle{empty}
	\openup 2.75\jot % double spaced
	\begin{center}
		\def\thanks##1{} % temporarily kill any thanks lines
		{\large \uppercase\expandafter{\@title} \par}
		{\large\lineskip .75em\begin{tabular}[t]{c}
		{\@author}\end{tabular}\par}
		{\large Under the supervision of Professor \@advisor \par}
		\vskip 1.5em
	\end{center}
}


% INTRODUCTION
%
%
% added by djj@mthcsc.wfu.edu on 3-27-90 as part of Wake Forest Thesis
% style
%
\newcommand{\introduction}{\pchapter{Introduction}}

% VITA
%
%
% added by djj@mthcsc.wfu.edu on 3-27-90 as part of Wake Forest Thesis
% style
%
%\newcommand{\vita}{\pchapter{Vita}}

% DEDICATION
%
%
% added djj@mthcsc.wfu.edu on 4-9-91 as part of Wake Forest Thesis
% style
%
\newcommand{\dedication}{\pchapter{Dedication}}

% ACKNOWLEDGEMENTS
%
%
\newcommand{\acknowledgements}{\pchapter{Acknowledgements}}
\newcommand{\acknowledgments}{\pchapter{Acknowledgments}}

% CHAPTERS
%
% pages following this command will be doublespaced with arabic numbering
\newcommand{\chapters}{
	\newpage
	\textspace
	\pagenumbering{arabic}
}

% APPENDIX
%
% The \appendix command must do the following:
%    -- reset the chapter counter to zero
%    -- set \@chapapp to Appendix (for messages)
%    -- redefine the chapter counter to produce appendix numbers
%    -- reset the section counter to zero
%    -- redefine the \chapter command if appendix titles and headings are
%       to look different from chapter titles and headings.

\newcommand{\appendix}{\par
  	\setcounter{chapter}{0}
  	\setcounter{section}{0}
	\renewcommand{\@chapapp}{Appendix}
  	\renewcommand{\thechapter}{\Alph{chapter}}}

%  VITA
%
%  djj -- May 17, 2012
%
%  The \vita coomand must do the following:
%    --  reset the chapter counter to zero
%    --  set \@chapapp to Curriculum Vitae
%    --  reset the section counter to zero
%    --  set \thechapter to the null string

\newcommand{\vita}{\par
        \setcounter{chapter}{0}
        \setcounter{section}{0}
        \renewcommand{\@chapapp}{Curriculum Vitae}
        \renewcommand{\thechapter}{}}


%    ****************************************
%    *                LISTS                 *
%    ****************************************
%

% The following commands are used to set the default values for the list
% environment's parameters. See the LaTeX manual for an explanation of
% the meanings of the parameters.  Defaults for the list environment are
% set as follows.  First, \rightmargin, \listparindent and \itemindent
% are set to 0pt.  Then, for a Kth level list, the command \@listK is
% called, where 'K' denotes 'i', 'ii', ... , 'vi'.  (I.e., \@listiii is
% called for a third-level list.)  By convention, \@listK should set
% \leftmargin to \leftmarginK.
%
% For efficiency, level-one list's values are defined at top level, and
% \@listi is defined to set only \leftmargin.

\leftmargini 2.5em
\leftmarginii 2.2em     % > \labelsep + width of '(m)'
\leftmarginiii 1.87em   % > \labelsep + width of 'vii.'
\leftmarginiv 1.7em     % > \labelsep + width of 'M.'
\leftmarginv 1em
\leftmarginvi 1em

\leftmargin\leftmargini
\labelwidth\leftmargini\advance\labelwidth-\labelsep
\labelsep .5em
\parsep 5pt plus 2.5pt minus 1pt

\newcommand{\@listi}{\leftmargin\leftmargini}

\newcommand{\@listii}{\leftmargin\leftmarginii
   	\labelwidth\leftmarginii\advance\labelwidth-\labelsep
   	\topsep 5pt plus 2.5pt minus 1pt
   	\parsep 2.5pt plus 1pt minus 1pt
   	\itemsep \parsep}

\newcommand{\@listiii}{\leftmargin\leftmarginiii
    	\labelwidth\leftmarginiii\advance\labelwidth-\labelsep
    	\topsep 2.5pt plus 1pt minus 1pt
    	\parsep \z@ \partopsep 1pt plus 0pt minus 1pt
    	\itemsep \topsep}

\newcommand{\@listiv}{\leftmargin\leftmarginiv
     	\labelwidth\leftmarginiv\advance\labelwidth-\labelsep}

\newcommand{\@listv}{\leftmargin\leftmarginv
     	\labelwidth\leftmarginv\advance\labelwidth-\labelsep}

\newcommand{\@listvi}{\leftmargin\leftmarginvi
     	\labelwidth\leftmarginvi\advance\labelwidth-\labelsep}


%    ****************************************
%    *                LISTS                 *
%    ****************************************
%

% ENUMERATE
%  Enumeration is done with four counters: enumi, enumii, enumiii
%  and enumiv, where enumN controls the numbering of the Nth level
%  enumeration.  The label is generated by the commands \labelenumi
%  ... \labelenumiv.  The expansion of \p@enumN\theenumN defines the
%  output of a \ref command.

\newcommand{\labelenumi}{\arabic{enumi}.}
\renewcommand{\theenumi}{\arabic{enumi}}

\newcommand{\labelenumii}{(\alph{enumii})}
\renewcommand{\theenumii}{\alph{enumii}}
\renewcommand{\p@enumii}{\theenumi}

\newcommand{\labelenumiii}{\roman{enumiii}.}
\renewcommand{\theenumiii}{\roman{enumiii}}
\renewcommand{\p@enumiii}{\theenumi(\theenumii)}

\newcommand{\labelenumiv}{\Alph{enumiv}.}
\renewcommand{\theenumiv}{\Alph{enumiv}}
\renewcommand{\p@enumiv}{\p@enumiii\theenumiii}

% ITEMIZE
% Itemization is controlled by four commands: \labelitemi, \labelitemii,
% \labelitemiii, and \labelitemiv, which define the labels of the various
% itemization levels.

\newcommand{\labelitemi}{$\bullet$}
\newcommand{\labelitemii}{\bf --}
\newcommand{\labelitemiii}{$\ast$}
\newcommand{\labelitemiv}{$\cdot$}


% VERSE
%   The verse environment is defined by making clever use of the
%   list environment's parameters.  The user types \\ to end a line.
%   This is implemented by \let'in \\ equal \@centercr.
%
\newcommand{\verse}{\let\\=\@centercr
  	\list{}{\itemsep\z@ \itemindent -1.5em\listparindent \itemindent
        \rightmargin\leftmargin\advance\leftmargin 1.5em}\item[]}
\let\endverse=\endlist

% QUOTATION
%   Fills lines
%   Indents paragraph
%
\newcommand{\quotation}{\list{}{\listparindent 1.5em
    	\itemindent\listparindent
    	\rightmargin\leftmargin\parsep 0pt plus 1pt}\item[]}
\let\endquotation=\endlist

% QUOTE -- same as quotation except no paragraph indentation,
%
\newcommand{\quote}{\list{}{\rightmargin\leftmargin}\item[]}
\let\endquote=\endlist

% DESCRIPTION
%
%  To change the formatting of the label, you must redefine
%  \descriptionlabel.

\newcommand{\descriptionlabel}[1]{\hspace\labelsep \bf #1}
\newcommand{\description}{\list{}{\labelwidth\z@ \itemindent-\leftmargin
       	\let\makelabel=\descriptionlabel}}

\let\enddescription=\endlist

\newdimen\descriptionmargin
\descriptionmargin=3em

%    ****************************************
%    *         OTHER ENVIRONMENTS           *
%    ****************************************
%
%
% THEOREM
% \@begintheorem ... \@endtheorem are the commands executed at the
% beginning and end of a (user-defined) theorem-like environment.
%
\renewcommand{\@begintheorem}[2]{\it \trivlist \item[\hskip 
	\labelsep{\bf #1\ #2}]}

\renewcommand{\@endtheorem}{\endtrivlist}

% EQUATION and EQNARRAY -- must follow \chapter definition
%
% \newcounter{equation}
%
\@addtoreset{equation}{chapter}   % Makes \chapter reset 'equation' counter.

\renewcommand{\theequation}{\thechapter.\arabic{equation}}

% \jot = 3pt      % Extra space added between lines of an eqnarray environment

% The macro \@eqnnum defines how equation numbers are to appear in equations.
%
% \newcommand{\@eqnnum}{(\theequation)}
%


% TITLEPAGE
%  In the normal environments, the titlepage environment does nothing but
%  start and end a page, and inhibit page numbers.  It also resets the
%  page number to two.

\newcommand{\titlepage}{\newpage\thispagestyle{empty}\c@page\tw@}

\def\endtitlepage{\newpage}
% defeat the setting to 0 of the page number for abstract pages.
\newcommand{\abstractpage}{\newpage}

\def\endabstractpage{\newpage}

% CODE
%  The code environment is essentially a 'quote'd verbatim environment.
%     \begin{code} ... \end{code}
%          ==>  {\quote\verbatim ... \endverbatim\endquote}
%  Page breaks are discouraged by \interlinepenalty=150.
%
\begingroup \catcode `|=0 \catcode `[=1
\catcode`]=2 \catcode`\{=12 \catcode `\}=12
\catcode`\\=12 |gdef|@xcode#1\end{code}[#1|end[code]]
|endgroup
\newcommand{\code}{\par\list{}{\parskip 0pt plus 1pt\normalbaselines
	\interlinepenalty=150}\item[]\@verbatim
	\frenchspacing\@vobeyspaces\@xcode}
\def\endcode{\endtrivlist\endlist}

% PROC
%  The proc environment is a 'quoted' tabbing environment with tabs set
%  at fixed distances
%  Page breaks are discouraged by \interlinepenalty=150.
%
\newcommand{\proctabs}{\qquad\=\qquad\=\qquad\=\qquad\=\qquad\
	=\qquad\=\qquad\=\kill}
\newenvironment{proc}{\quote\interlinepenalty=150\relax\tt
	\tabbing\proctabs}{\endtabbing\endquote\leavevmode}

% ARRAY AND TABULAR
%

\arraycolsep 5pt     % Half the space between columns in an array environment.
\tabcolsep 6pt       % Half the space between columns in a tabular environment.
\arrayrulewidth .4pt % Width of rules in array and tabular environment.
\doublerulesep 2pt   % Space between adjacent rules in array or tabular env.

% TABBING
%
\tabbingsep \labelsep   % Space used by the \' command.  (See LaTeX manual.)

% MINIPAGE
%  \@minipagerestore is called upon entry to a minipage environment to
%  set up things that are to be handled differently inside a minipage
%  environment. In the current styles, it does nothing.
%
% \skip\@mpfootins : plays same role for footnotes in a minipage as
%                    \skip\footins does for ordinary footnotes

\skip\@mpfootins = \skip\footins

% FRAMEBOX
%
\fboxsep = 3pt    % Space left between box and text by \fbox and \framebox.
\fboxrule = .4pt  % Width of rules in box made by \fbox and \framebox.

%    ****************************************
%    *         TABLE OF CONTENTS, ETC.      *
%    ****************************************
%
% A \subsection command writes a
%       \contentsline{subsection}{TITLE}{PAGE}
% command on the .toc file, where TITLE contains the contents of the
% entry and PAGE is the page number.  If subsections are being numbered,
% then TITLE will be of the form
%       \numberline{NUM}{HEADING}
% where NUM is the number produced by \thesubsection.  Other sectioning
% commands work similarly.
%
% A \caption command in a 'figure' environment writes
%    \contentsline{figure}{\numberline{NUM}{CAPTION}}{PAGE}
% on the .lof file, where NUM is the number produced by \thefigure and
% CAPTION is the figure caption.  It works similarly for a 'table' environment.
%
% The command \contentsline{NAME} expands to \l@NAME.  So, to specify
% the table of contents, we must define \l@chapter, \l@section,
% \l@subsection, ... ; to specify the list of figures, we must define
% \l@figure; and so on.  Most of these can be defined with the
% \@dottedtocline command, which works as follows.
%
% \@dottedtocline{LEVEL}{INDENT}{NUMWIDTH}{TITLE}{PAGE}
%    LEVEL    : An entry is produced only if LEVEL < or = value of
%               'tocdepth' counter.  Note, \chapter is level 0, \section
%               is level 1, etc.
%    INDENT   : The indentation from the outer left margin of the start of
%               the contents line.
%    NUMWIDTH : The width of a box in which the section number is to go,
%               if TITLE includes a \numberline command.
%
% This command uses the following three parameters, which are set
% with a \newcommand (so em's can be used to make them depend upon the font).
%   \@pnumwidth : The width of a box in which the page number is put.
%   \@tocrmarg  : The right margin for multiple line entries.  One
%                 wants \@tocrmarg > or = \@pnumwidth
%   \@dotsep    : Separation between dots, in mu units.  Should be \def'd to
%                 a number like 2 or 1.7

\newcommand{\@pnumwidth}{1.55em}
\newcommand{\@tocrmarg} {2.55em}
\newcommand{\@dotsep}{4.5}
\setcounter{tocdepth}{2}


% TABLEOFCONTENTS
%

\newcommand{\tableofcontents}{
  	\normalbaselines
 	\chapter*{Table of Contents}
  	\thispagestyle{headings}
  	\typeout{Table of Contents}
  	\markboth{CONTENTS}{CONTENTS}
  	\@starttoc{toc}}

% NOTE:  unbolded to comply with WFU Thesis style (May, 2012)

\newcommand{\l@chapter}[2]{\pagebreak[3]
	\vskip 1.0em plus 1pt        % space above chapter line
   	\@tempdima 1.5em       	     % width of box holding chapter number
   	\begingroup
     		\parindent \z@ \rightskip \@pnumwidth
     		\parfillskip -\@pnumwidth
%     		\bf                  % Boldface.
     		\leavevmode          % TeX command to enter horizontal mode.
     		#1\dotfill \hbox to\@pnumwidth{\hss #2}\par
   	\endgroup}

\newcommand{\l@section}{\@dottedtocline{1}{1.5em}{2.3em}}
\newcommand{\l@subsection}{\@dottedtocline{2}{3.8em}{3.2em}}
\newcommand{\l@subsubsection}{\@dottedtocline{3}{7.0em}{4.1em}}
\newcommand{\l@paragraph}{\@dottedtocline{4}{10em}{5em}}
\newcommand{\l@subparagraph}{\@dottedtocline{5}{12em}{6em}}

% LIST OF FIGURES
%
\newcommand{\listoffigures}{
	\normalbaselines
	\pchapter{List of Figures}
	\markboth{LIST OF FIGURES}{LIST OF FIGURES}
	\@starttoc{lof}
}

\newcommand{\l@figure}{\@dottedtocline{1}{1.5em}{2.3em}}

% LIST OF TABLES
%
\newcommand{\listoftables}{
	\normalbaselines
	\pchapter{List of Tables}
	\markboth{LIST OF TABLES}{LIST OF TABLES}
	\@starttoc{lot}
}

\let\l@table\l@figure


% GREMLIN
%
%  \Gremlin{includefile}{width}{height}{F/E}
%     Warning: ignoring the Flyback/End option
%
\newcommand{\Gremlin}[4]{
	\leavevmode\hbox to #2{\vbox to #3{\special{#1}\vfil}\hfil}
}

%    ****************************************
%    *             BIBLIOGRAPHY             *
%    ****************************************
%
%  TEMPORARY DEFINITIONS:

\newcommand{\thebibliography}[1]{
  	\chapter*{Bibliography\markboth{BIBLIOGRAPHY}{BIBLIOGRAPHY}}
  	\list{[\arabic{enumi}]}{\settowidth\labelwidth{[#1]}
	\leftmargin\labelwidth
    	\advance\leftmargin\labelsep\usecounter{enumi}}
  	\addcontentsline{toc}{chapter}{Bibliography}
  	\newcommand{\newblock}{\hskip .11em plus .33em minus -.07em}
  	\sloppy
  	\sfcode`\.=1000\relax
}


\let\endthebibliography=\endlist

% \newcommand{\@biblabel}[1]{[#1]\hfill}  % Produces the label for a 
					  % \bibitem[...] command.
% \newcommand{\@cite}[1]{[#1]}            % Produces the output of the
					  % \cite command.


%    ****************************************
%    *              THE INDEX               *
%    ****************************************
%
% THE INDEX ENVIRONMENT
% Produces double column format, with each paragraph a separate entry.
% The user commands \item, \subitem and \subsubitem are used to
% produce the entries, and \indexspace adds an extra vertical space
% that's the right size to put above the first entry with a new letter
% of the alphabet.

\newcommand{\theindex}{
	\columnseprule \z@
	\columnsep 35pt\twocolumn[\@makeschapterhead{Index}]
    	\markboth{INDEX}{INDEX}\thispagestyle{plain}\parindent\z@
    	\parskip\z@ plus .3pt\relax\let\item\@idxitem}

\newcommand{\@idxitem}{\par\hangindent 40pt}

\newcommand{\subitem}{\par\hangindent 40pt \hspace*{20pt}}

\newcommand{\subsubitem}{\par\hangindent 40pt \hspace*{30pt}}

\def\endtheindex{\onecolumn}

\newcommand{\indexspace}{\par \vskip 10pt plus 5pt minus 3pt\relax}

%    ****************************************
%    *             FOOTNOTES                *
%    ****************************************
%
% \footnoterule is a macro to draw the rule separating the footnotes from
% the text.  It should take zero vertical space, so it needs a negative
% skip to compensate for any positive space taken by the rule.  (See
% PLAIN.TEX.)

\renewcommand{\footnoterule}{\kern-1\p@
  	\hrule width .4\columnwidth
  	\kern .6\p@}                 % The \hrule has default height of .4pt .

% \newcounter{footnote}
\@addtoreset{footnote}{chapter}  % Numbers footnotes within chapters

%   \@makefntext{NOTE} :
%        Must produce the actual footnote, using \@thefnmark as the mark
%        of the footnote and NOTE as the text.  It is called when effectively
%        inside a \parbox of width \columnwidth (i.e., with \hsize =
%        \columnwidth).  The following macro indents
%        all lines of the footnote by 10pt, and indents the first line of
%        a new paragraph by 1em.  To change these dimensions, just substitute
%        the desired value for '10pt' [in both places] or '1em'.
%        The mark is flushright against the footnote.
%
%        It is much easier to write a macro in which the footnote text is
%        set like an ordinary text paragraph, with no indentation except
%        on the first line of a paragraph, and the first line of the
%        footnote.  In that case, all the macro must do is set \parindent
%        to the appropriate value for succeeding paragraphs and put the
%        proper indentation before mark.

\long\def\@makefntext#1{\@setpar{\@@par\@tempdima \hsize
  	\advance\@tempdima-10pt\parshape \@ne 10pt \@tempdima}\par
  	\parindent 1em\noindent \hbox to \z@{\hss$^{\@thefnmark}$}#1}


% \@makefnmark : A macro to generate the footnote marker that goes
%    in the text.  Default used.
%


%    ****************************************
%    *         FIGURES AND TABLES           *
%    ****************************************
%
% Float placement parameters.  See LaTeX manual for their definition.
%
\setcounter{topnumber}{2}
\renewcommand{\topfraction}{.7}
\setcounter{bottomnumber}{1}
\renewcommand{\bottomfraction}{.3}
\setcounter{totalnumber}{3}
\renewcommand{\textfraction}{.2}
\renewcommand{\floatpagefraction}{.5}
\setcounter{dbltopnumber}{2}
\renewcommand{\dbltopfraction}{.7}
\renewcommand{\dblfloatpagefraction}{.5}

%
% fix a BUG? in latex.tex
%
\let\@@addtocurcol\@addtocurcol
\renewcommand{\@addtocurcol}{{\@floatplacement\@@addtocurcol}}

% \@makecaption{NUMBER}{TEXT} : Macro to make a figure or table caption.
%      NUMBER : Figure or table number--e.g., 'Figure 3.2'
%      TEXT   : The caption text.
%  Macro should be called inside a \parbox of right width, with \normalsize.

\long\def\@makecaption#1#2{
   	\vskip 10pt
   	\setbox\@tempboxa\hbox{#1: #2}
   	\ifdim \wd\@tempboxa >\hsize   % IF longer than one line:
            \unhbox\@tempboxa\par      %   THEN set as ordinary paragraph.
     	\else                          %   ELSE  center.
            \hbox to\hsize{\hfil\box\@tempboxa\hfil}
   	\fi}

% To define a float of type TYPE (e.g., TYPE = figure), the document style
% must define the following.
%
%  \fps@TYPE   : The default placement specifier for floats of type TYPE.
%
%  \ftype@TYPE : The type number for floats of type TYPE.  Each TYPE has
%                associated a unique positive TYPE NUMBER, which is a power
%                of two.  E.g., figures might have type number 1, tables type
%                number 2, programs type number 4, etc.
%
%  \ext@TYPE   : The file extension indicating the file on which the
%                contents list for float type TYPE is stored.  For example,
%                \ext@figure = 'lof'.
%
%  \fnum@TYPE  : A macro to generate the figure number for a caption.
%                For example, \fnum@TYPE == Figure \thefigure.
%
%  The actual float-making environment commands--e.g., the commands
%  \figure and \endfigure--are defined in terms of the macros \@float
%  and \end@float, which are described below.
%
%  \@float{TYPE}[PLACEMENT] : Macro to begin a float environment for a
%     single-column float of type TYPE with PLACEMENT as the placement
%     specifier.  The default value of PLACEMENT is defined by \fps@TYPE.
%     The environment is ended by \end@float.
%     E.g., \figure == \@float{figure}, \endfigure == \end@float.


% FIGURE

\newcounter{figure}[chapter]
\renewcommand{\thefigure}{\thechapter.\@arabic\c@figure}

\newcommand{\fps@figure}{tbp}
\newcommand{\ftype@figure}{1}
\newcommand{\ext@figure}{lof}
\newcommand{\fnum@figure}{Figure \thefigure}
\newcommand{\figure}{\@float{figure}}
\let\endfigure\end@float
\@namedef{figure*}{\@dblfloat{figure}}
\@namedef{endfigure*}{\end@dblfloat}

% TABLE
%
\newcounter{table}[chapter]
\renewcommand{\thetable}{\thechapter.\@arabic\c@table}

\newcommand{\fps@table}{tbp}
\newcommand{\ftype@table}{2}
\newcommand{\ext@table}{lot}
\newcommand{\fnum@table}{Table \thetable}
\newcommand{\table}{\@float{table}}
\let\endtable\end@float
\@namedef{table*}{\@dblfloat{table}}
\@namedef{endtable*}{\end@dblfloat}

%    ****************************************
%    *         TITLE AND ABSTRACT           *
%    ****************************************
%


% \department{department name} -- user macro to declare name of department.
\newcommand{\department}[1]{\gdef\@department{#1}}

% \advisor{advisor name} -- user macro to declare name of advisor/supervisor.
\newcommand{\advisor}[1]{\gdef\@advisor{#1}}

% \chairperson{chairperson name} -- user macro to declare name of committee
%                                   chairperson
\newcommand{\chairperson}[1]{\gdef\@chairperson{#1}}

% Set default values for department name, advisor name, and current year.
\department{??DEPARTMENT??}
\advisor{??ADVISOR??}
\chairperson{??CHAIRPERSON??}
\date{\number\year}

% The macros and counters between here and the definition of \maketitle help 
% \maketitle handle the signature lines.

% \@member{member name} -- internal macro to print a line to sign on, and the
%	name of a committee member (sole argument) below the line.
%  signing line removed (2011)
\newcommand{\@member}[1]{
%	#1  \hfill {\rule{2.375in}{.5pt}} \par
        #1  \par
}

\newcounter{@members}		% Number of member macros defined so far.
\newcommand{\@membera}{\@member{\@advisor, Advisor}}
\newcommand{\@memberb}{\@member{\@chairperson, Chair}}
\setcounter{@members}{2}	% Advisor is always 1st member and chairperson
                                % is always second member

% @minmembers -- holds mimumum number of member macros needed.  It is set in
% \maketitle.
\newcounter{@minmembers}

% \member{member name} -- user macro to declare the name of a committee member.
%	May be used up to 7 times.  Should not be used for advisor.
\newcommand{\member}[1]{
  	\addtocounter{@members}{1}	% Inc. num of members to 
					% include this new one.
  	\ifcase\value{@members} 
    	    \typeout{Thesis format internal error.  Too few committee members.} 
 \or
    	    \typeout{Thesis format internal error.  Advisor left off committee.}
 \or
    	    \typeout{Thesis format internal error.  Chairperson left off committee.} \or
    	    \newcommand{\@memberc}{\@member{#1}} \or
    	    \newcommand{\@memberd}{\@member{#1}} \or
    	    \newcommand{\@membere}{\@member{#1}} \or
    	    \newcommand{\@memberf}{\@member{#1}} \or
    	    \newcommand{\@memberg}{\@member{#1}} \or
    	    \newcommand{\@memberh}{\@member{#1}} 
  	\else
    	    \typeout{Error using thesis format.  Too many committee members.}
  	\fi
}


%- Define an abstract macro to provide an abstract page that exactly
%- duplicates the title page, EXCEPT for the words ``an abstract'' vice
%- ``a dissertation.''
%-
\newcommand{\maketitle}{		% Define macro to print title page.
  	\typeout{Title page}
  	\begin{titlepage}
    	  	\let\footnotesize\small % Local definitions to make 
					% \thanks produce
    	  	\let\footnoterule\relax % \normalsize footnotes with
					%  no separating rule.
    		\setcounter{page}{1}
    	\null
    	\parskip=0pt
    	\vspace*{60pt}
    	\begin{center}
      		\vglue -.25in   % this makes top margin 2in
      		{\normalsize\uppercase\expandafter{\@title} \par}
      		\vglue 5ex
      		{\normalsize BY \par}
      		\vglue 5ex
      		{\normalsize \uppercase\expandafter{\@author}\par\par}
      		\vglue 3ex
    	\end{center}


    	\begin{center}
		\normalsize 
		\ifPhD
	  	A Dissertation submitted to the Graduate Faculty of\par
		\else
	  	A Thesis Submitted to the Graduate Faculty of \par
		\fi
	  	\vglue 2ex
	  	\mbox{WAKE$\:$FOREST$\:$UNIVERSITY$\:$GRADUATE$\:$SCHOOL$\:$OF$\:$ARTS$\:$AND$\:$SCIENCES} \par
	  	\vglue 2ex
	  	in Partial Fulfillment of the Requirements \par
	  	\vglue 2ex
	  	for the Degree of \par
	  	\vglue 2ex
	  	\ifMS
	    	MASTER OF SCIENCE \par
		\fi
	  	\ifAM
	    	MASTER OF ARTS \par
		\fi
		\ifPhD
		Doctor of Philosophy \par
	  	\fi
		\vglue 2ex
		\@department \par
		\vglue 2ex
		\@date \par
		\vglue 2ex
		Winston-Salem, North Carolina
    	\end{center}
    
    
 %   	\vglue 3ex plus 1fill  % fill to keep bottom of page from creeping up
	\vglue 3ex
    	\ifPhD
      	\setcounter{@minmembers}{5} % PhD committees have at least 5 members
    	\else
      	\setcounter{@minmembers}{3} % Masters committees have at least 3 members
    	\fi
%    	\parbox{5.5in} {
	\begin{center}
      		Approved By:\par
      		\vglue 2ex
      		\@membera
      		\vglue 2ex
      		%Examining Committee: \par
      		%\vglue 2ex	
 %     		{ 
	% @membersprinted counts the number of committee members printed so far
	\newcounter{@membersprinted}\setcounter{@membersprinted}{1}
	% Print the committee members defined by the user.
	\loop
	\ifnum \value{@membersprinted} < \value{@members}
	  \addtocounter{@membersprinted}{1}
	  \ifcase\value{@membersprinted}
	    \typeout{Thesis style--internal error printing members.} \or
	    \@membera \or
	    \@memberb \or
	    \@memberc \or
	    \@memberd \or
	    \@membere \or
	    \@memberf \or
	    \@memberg \or
	    \@memberh
	  \else \typeout{Thesis style--internal error printing members.}
	  \fi
	  \vglue 1.5ex
	\repeat
	% If the minimum number of signature lines haven't been printed,
	% print some more.
	\loop
	\ifnum \value{@membersprinted} < \value{@minmembers}
	  \addtocounter{@membersprinted}{1}
	  \@member{}
	\repeat
 %     \end{center}
 %   }\hfill
    \vglue 2ex
    \end{center}
 
  \end{titlepage}
  \setcounter{footnote}{0} % Footnotes start at zero again.
  \let\thanks\relax
%%% SAVE THESE \gdef\@thanks{}\gdef\@department{}\gdef\@author{}\gdef\@title{}
  \let\maketitle\relax % only ONCE!?
}
%% Makeabstract command
\newcommand{\makeabstract}{		% Define macro to print title page.
  	\typeout{Abstract page}
  	\begin{abstractpage}
    	    \let\footnotesize\small % Local definitions to make \thanks produce
    	    \let\footnoterule\relax % \normalsize footnotes with no separating rule.
    	    \setcounter{page}{0}
    	    \null
    	    \parskip=0pt
    	    \begin{center}
      		\vglue -.25in   % this makes top margin 2in
      		{\Large \uppercase\expandafter{\@title} \par}
      		\vglue 2ex
      		{\normalsize by \par}
      		\vglue 2ex
      		{\large\@author\par\par}
      		\vglue 3ex
      		{\normalsize Department of \@department\par}
      		{\normalsize Wake Forest University\par}
    	\end{center}

    \vglue 3ex plus 1fill  % fill to keep bottom of page from creeping up
    \ifPhD
      \setcounter{@minmembers}{5} % PhD committees have at least 5 members
    \else
      \setcounter{@minmembers}{3} % Masters committees have at least 3 members
    \fi
    \parbox{2.8in}{
      Date:\hfill\rule[-.4ex]{2.35in}{.5pt}\par
      Approved:\par
      { \centering\normalsize
	% @membersprinted counts the number of committee members printed so far
	\setcounter{@membersprinted}{0}
	% Print the committee members defined by the user.
	\loop
	\ifnum \value{@membersprinted} < \value{@members} 
	  \addtocounter{@membersprinted}{1}
	  \ifcase\value{@membersprinted}
	    \typeout{Thesis style--internal error printing members.} \or
	    \@membera \or
	    \@memberb \or
	    \@memberc \or
	    \@memberd \or
	    \@membere \or
	    \@memberf \or
	    \@memberg \or
	    \@memberh
	  \else \typeout{Thesis style--internal error printing members.}
	  \fi
	\repeat
	% If the minimum number of signature lines haven't been printed,
	% print some more.
	\loop
	\ifnum \value{@membersprinted} < \value{@minmembers}
	  \addtocounter{@membersprinted}{1}
	  \@member{}
	\repeat
      }
    }\hfill
    \vglue 2ex

    \begin{center}
      {\normalsize 
	An abstract of a dissertation submitted in partial\par
        fulfillment of the requirements for the degree \par
        of Doctor of Philosophy in the Department of \par
        \@department\ in the Graduate School of \par
	Wake Forest University \par
	\vglue 0.75em
	\@date \par
      }
    \end{center}
    \par\vglue 1ex\@thanks\null
  \end{abstractpage}
  \setcounter{footnote}{0} % Footnotes start at zero again.
  \let\thanks\relax
%%% SAVE THESE \gdef\@thanks{}\gdef\@department{}\gdef\@author{}\gdef\@title{}
  \let\makeabstract\relax % only ONCE!?
}
%% End of the make abstract

% Old junk; probably should be deleted.
%\newcommand{\Type}[1]{\gdef\@Type{#1}}
%\Type{Dissertation}

%    ****************************************
%    *            PAGE STYLES               *
%    ****************************************
%
% The page style 'foo' is defined by defining the command \ps@foo.  This
% command should make only local definitions.  There should be no stray
% spaces in the definition, since they could lead to mysterious extra
% spaces in the output.
%
% The \ps@... command defines the macros \@oddhead, \@oddfoot,
% \@evenhead, and \@evenfoot to define the running heads and
% feet---e.g., \@oddhead is the macro to produce the contents of the
% heading box for odd-numbered pages.  It is called inside an \hbox of
% width \textwidth.
%
% To make headings determined by the sectioning commands, the page style
% defines the commands \chaptermark, \sectionmark, ... , where
% \chaptermark{TEXT} is called by \chapter to set a mark, and so on.
% The \...mark commands and the \...head macros are defined with the
% help of the following macros.  (All the \...mark commands should be
% initialized to no-ops.)
%
% MARKING CONVENTIONS:
% LaTeX extends TeX's \mark facility by producing two kinds of marks
% a 'left' and a 'right' mark, using the following commands:
%     \markboth{LEFT}{RIGHT} : Adds both marks.
%     \markright{RIGHT}      : Adds a 'right' mark.
%     \leftmark  : Used in the \@oddhead, \@oddfoot, \@evenhead or \@evenfoot
%                  macro, gets the current 'left'  mark.  Works like TeX's
%                  \botmark command.
%     \rightmark : Used in the \@oddhead, \@oddfoot, \@evenhead or \@evenfoot
%                  macro, gets the current 'right'  mark.  Works like TeX's
%                  \firstmark command.
% The marking commands work reasonably well for right marks 'numbered
% within' left marks--e.g., the left mark is changed by a \chapter command and
% the right mark is changed by a \section command.  However, it does
% produce somewhat anomalous results if two \bothmark's occur on the same page.
%

\mark{{}{}}   % Initializes TeX's marks

\renewcommand{\ps@empty}{\renewcommand\@oddhead{\versiontag\hfil} 
	\renewcommand{\@oddfoot}{}
	\renewcommand{\@evenhead}{\versiontag\hfil}
	\renewcommand{\@evenfoot}{}}

%
% Definiton of 'plain' page style.
%   This page style is used on 2nd and subsequent pages of a chapter

\renewcommand{\ps@plain}{%
%	\renewcommand{\@oddhead}{\versiontag\hfil\rm\thepage}
	\renewcommand{\@oddhead}{\versiontag\hfil}
	\renewcommand{\@oddfoot}{\hfil\rm\thepage\hfil}
	\let\@evenhead\@oddhead
	\let\@evenfoot\@oddfoot}

% Definition of 'headings' page style
%  Note the use of ##1 for parameter of \def\chaptermark inside the
%  \newcommand{\ps@headings}.
%

\newcommand{\ps@headings}{%
	\renewcommand\@oddfoot{\hfil\rm\thepage\hfil}%  centered page number
	\let\@evenfoot\@oddfoot
	\renewcommand\@oddhead{\versiontag\hfil}
	\let\@evenhead\@oddhead
	\def\chaptermark##1{\markright {\uppercase{\ifnum 
		\c@secnumdepth >\m@ne \@chapapp\ \thechapter. \ \fi ##1}}}}

% Definition of 'myheadings' page style.  The user's \head command
% changes the definition of \@lhead and/or \@rhead.
%
\newcommand{\ps@myheadings}{
	\renewcommand{\@oddhead}{\hbox{}\sl\markright \hfil \rm\thepage}%
	\renewcommand{\@oddfoot}{}
	\renewcommand{\@evenhead}{\rm \thepage\hfil\sl\markleft \hbox {}}%
	\renewcommand{\@evenfoot}{}
	\renewcommand{\chaptermark}[1]{}
	\renewcommand{\sectionmark}[1]{}}


%    ****************************************
%    *            MISCELLANEOUS             *
%    ****************************************
%
% DATE
%
\newcommand{\today}{\ifcase\month\or
  	January\or February\or March\or April\or May\or June\or
  	July\or August\or September\or October\or November\or December\fi
  	\space\number\day, \number\year}

%    ****************************************
%    *           INITIALIZATION             *
%    ****************************************
%
% Default initializations

\ps@plain                   % 'plain' page style
\pagenumbering{roman}       % roman page numbers
\onecolumn                  % Single-column.
\raggedbottom               % Ragged bottom

%    ****************************************
%    *            MODIFIED 4/9/2208         *
%    ****************************************
% Changed the spacing between the title, by, and authors name to tripple spacing
% Added code to force the authors name to all upper caps

%    ****************************************
%    *            USAGE NOTES 4/9/2208      *
%    ****************************************
% Changing between single and double spacing is an issue if an appendix is used as
% the bibliography must be single spaced and the appendix, located after the bib, must be double spaced.
% \normalbaselines forces single space
% \textspace forces double space
%
% the typical included packages are:
%\usepackage{graphicx} 		%pretty pictures
%\usepackage{float}			%counter used to force figures in place
%\usepackage{hyperref}		%click on reference to auto-scroll to table/figure
%\usepackage{subfigure}  	%used to place figures side-by-side
