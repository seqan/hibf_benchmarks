"""Functions for creating CSS and HTML components for the Bokeh plot."""

import os
from bokeh.models import Div


def create_vercel_div():
    """Creates a div containing the Vercel logo."""
    vercel_logo = """
        <svg width="157" height="30" viewBox="0 0 209 40" fill="none" xmlns="http://www.w3.org/2000/svg">
        <path d="M0 5C0 2.23858 2.23858 0 5 0H204C206.761 0 209 2.23858 209 5V35C209 37.7614 206.761 40 204 40H5C2.23858 40 0 37.7614 0 35V5Z" fill="black"/>
        <path fill-rule="evenodd" clip-rule="evenodd" d="M20 13L28 27H12L20 13Z" fill="white"/>
        <line x1="40.5" y1="2.18556e-08" x2="40.5" y2="40" stroke="#333333"/>
        <path d="M53.2784 26H55.0341V21.9091H57.4205C60.1193 21.9091 61.4545 20.2784 61.4545 18.1307C61.4545 15.9886 60.1307 14.3636 57.4261 14.3636H53.2784V26ZM55.0341 20.4205V15.8693H57.2386C58.9773 15.8693 59.6875 16.8125 59.6875 18.1307C59.6875 19.4489 58.9773 20.4205 57.2614 20.4205H55.0341ZM66.9432 26.1761C69.4034 26.1761 71.0114 24.375 71.0114 21.6761C71.0114 18.9602 69.4034 17.1591 66.9432 17.1591C64.483 17.1591 62.875 18.9602 62.875 21.6761C62.875 24.375 64.483 26.1761 66.9432 26.1761ZM66.9489 24.75C65.3409 24.75 64.5909 23.3466 64.5909 21.6705C64.5909 20 65.3409 18.5795 66.9489 18.5795C68.5455 18.5795 69.2955 20 69.2955 21.6705C69.2955 23.3466 68.5455 24.75 66.9489 24.75ZM74.5341 26H76.2614L78.0341 19.6989H78.1648L79.9375 26H81.6705L84.233 17.2727H82.4773L80.7784 23.6534H80.6932L78.9886 17.2727H77.233L75.517 23.6818H75.4318L73.7216 17.2727H71.9659L74.5341 26ZM89.3409 26.1761C91.2443 26.1761 92.5909 25.2386 92.9773 23.8182L91.3693 23.5284C91.0625 24.3523 90.3239 24.7727 89.358 24.7727C87.9034 24.7727 86.9261 23.8295 86.8807 22.1477H93.0852V21.5455C93.0852 18.392 91.1989 17.1591 89.2216 17.1591C86.7898 17.1591 85.1875 19.0114 85.1875 21.6932C85.1875 24.4034 86.767 26.1761 89.3409 26.1761ZM86.8864 20.875C86.9545 19.6364 87.8523 18.5625 89.233 18.5625C90.5511 18.5625 91.4148 19.5398 91.4205 20.875H86.8864ZM94.9702 26H96.669V20.6705C96.669 19.5284 97.5497 18.7045 98.7543 18.7045C99.1065 18.7045 99.5043 18.767 99.6406 18.8068V17.1818C99.4702 17.1591 99.1349 17.142 98.919 17.142C97.8963 17.142 97.0213 17.7216 96.7031 18.6591H96.6122V17.2727H94.9702V26ZM104.56 26.1761C106.463 26.1761 107.81 25.2386 108.196 23.8182L106.588 23.5284C106.281 24.3523 105.543 24.7727 104.577 24.7727C103.122 24.7727 102.145 23.8295 102.099 22.1477H108.304V21.5455C108.304 18.392 106.418 17.1591 104.44 17.1591C102.009 17.1591 100.406 19.0114 100.406 21.6932C100.406 24.4034 101.986 26.1761 104.56 26.1761ZM102.105 20.875C102.173 19.6364 103.071 18.5625 104.452 18.5625C105.77 18.5625 106.634 19.5398 106.639 20.875H102.105ZM113.456 26.1705C115.047 26.1705 115.672 25.1989 115.979 24.642H116.121V26H117.78V14.3636H116.081V18.6875H115.979C115.672 18.1477 115.092 17.1591 113.467 17.1591C111.359 17.1591 109.808 18.8239 109.808 21.6534C109.808 24.4773 111.337 26.1705 113.456 26.1705ZM113.831 24.7216C112.314 24.7216 111.524 23.3864 111.524 21.6364C111.524 19.9034 112.297 18.6023 113.831 18.6023C115.314 18.6023 116.109 19.8125 116.109 21.6364C116.109 23.4716 115.297 24.7216 113.831 24.7216ZM124.575 26H126.234V24.642H126.376C126.683 25.1989 127.308 26.1705 128.899 26.1705C131.013 26.1705 132.547 24.4773 132.547 21.6534C132.547 18.8239 130.99 17.1591 128.882 17.1591C127.263 17.1591 126.678 18.1477 126.376 18.6875H126.274V14.3636H124.575V26ZM126.24 21.6364C126.24 19.8125 127.036 18.6023 128.518 18.6023C130.058 18.6023 130.831 19.9034 130.831 21.6364C130.831 23.3864 130.036 24.7216 128.518 24.7216C127.058 24.7216 126.24 23.4716 126.24 21.6364ZM135.216 29.25C136.619 29.25 137.511 28.517 138.011 27.1648L141.619 17.2898L139.784 17.2727L137.574 24.0455H137.483L135.273 17.2727H133.455L136.648 26.1136L136.438 26.6932C136.006 27.8239 135.398 27.9261 134.466 27.6705L134.057 29.0625C134.261 29.1591 134.705 29.25 135.216 29.25ZM149.426 14.3636H146.693L150.71 26H153.881L157.892 14.3636H155.165L152.347 23.2045H152.239L149.426 14.3636ZM162.224 26.1705C164.384 26.1705 165.838 25.1193 166.179 23.5L163.94 23.3523C163.696 24.017 163.071 24.3636 162.264 24.3636C161.054 24.3636 160.287 23.5625 160.287 22.2614V22.2557H166.23V21.5909C166.23 18.625 164.435 17.1591 162.128 17.1591C159.56 17.1591 157.895 18.983 157.895 21.6761C157.895 24.4432 159.537 26.1705 162.224 26.1705ZM160.287 20.7557C160.338 19.7614 161.094 18.9659 162.168 18.9659C163.219 18.9659 163.946 19.7159 163.952 20.7557H160.287ZM167.81 26H170.23V21.0625C170.23 19.9886 171.014 19.25 172.082 19.25C172.418 19.25 172.878 19.3068 173.105 19.3807V17.233C172.889 17.1818 172.588 17.1477 172.344 17.1477C171.366 17.1477 170.565 17.7159 170.247 18.7955H170.156V17.2727H167.81V26ZM177.893 26.1705C180.217 26.1705 181.678 24.8068 181.791 22.8011H179.507C179.365 23.733 178.751 24.2557 177.922 24.2557C176.791 24.2557 176.058 23.3068 176.058 21.6364C176.058 19.9886 176.797 19.0455 177.922 19.0455C178.808 19.0455 179.376 19.6307 179.507 20.5H181.791C181.689 18.483 180.161 17.1591 177.882 17.1591C175.234 17.1591 173.598 18.9943 173.598 21.6705C173.598 24.3239 175.206 26.1705 177.893 26.1705ZM187.318 26.1705C189.477 26.1705 190.932 25.1193 191.273 23.5L189.034 23.3523C188.79 24.017 188.165 24.3636 187.358 24.3636C186.148 24.3636 185.381 23.5625 185.381 22.2614V22.2557H191.324V21.5909C191.324 18.625 189.528 17.1591 187.222 17.1591C184.653 17.1591 182.989 18.983 182.989 21.6761C182.989 24.4432 184.631 26.1705 187.318 26.1705ZM185.381 20.7557C185.432 19.7614 186.188 18.9659 187.261 18.9659C188.312 18.9659 189.04 19.7159 189.045 20.7557H185.381ZM195.324 14.3636H192.903V26H195.324V14.3636Z" fill="white"/>
        </svg>
    """
    vercel_div = Div(
        text=f"<a href='https://vercel.com/?utm_source=seqan&utm_campaign=oss'><svg height=30px>{vercel_logo}</svg></a>",
        styles={
            "margin": "0px",
            "display": "block" if "VERCEL" in os.environ else "none",
        },
    )
    return vercel_div


def create_latex_text():
    """Creates a div containing the description of the plot."""
    return """
        <div>
            <strong>Parameters:</strong><br>
            <ul>
                <li><span>Alpha (\\(\\alpha\\)): Influences the ratio of merged bins and split bins</span></li>
                <li><span>Hash: The number of \\(hash\\) functions for Bloom Filters</span></li>
                <li><span>\\(k\\)-mer: Choosing window and k-mer size</span></li>
                <li><span>\\(fpr\\): Sets an upper bound for Bloom Filter false positives</span></li>
                <li><span>Estimate Union (\\(U\\)): Algorithm estimates sequence similarity between input data</span></li>
                <li><span>Rearrangement (\\(U+R\\)): Change order of sequences based on their estimated similarity</span></li>
            </ul>
            <strong>Default parameters:</strong><br>
            <ul>
                <li><span>\\(\\alpha = 1.2\\)</span></li>
                <li><span>\\(t_{\\text{max}} = 192\\)</span></li>
                <li><span>\\(\\text{hash} = 2\\)</span></li>
                <li><span>\\((w,k)=(32,32)\\)</span></li>
                <li><span>\\(\\text{r-fpr} = 0.5\\)</span></li>
                <li><span>\\(\\text{fpr} = 0.05\\)</span></li>
                <li><span>\\(U = True\\)</span></li>
                <li><span>\\(R = True\\)</span></li>
            </ul>
            <strong>Legend:</strong><br>
            <ul>
              <li><span>Determine query length: Time to determine the length of the query</span></li>
              <li><span>Queryfile IO: Time to read the query file</span></li>
              <li><span>Load index: Time to load the index</span></li>
              <li><span>Compute minimizer (avg): Average time to compute minimizer per thread</span></li>
              <li><span>Query IBF (avg): Average time to query the IBF per thread</span></li>
              <li><span>Generate results (avg): Average time to generate results per thread</span></li>
              <li><span>Level 0: Size of the first level of the index</span></li>
              <li><span>Level 1: Size of the second level of the index</span></li>
              <li><span>Level 2: Size of the third level of the index</span></li>
              <li><span>Level 3: Size of the fourth level of the index</span></li>
              <li><span>Avg load factor: Average load factor of the index</span></li>
            </ul>
        </div>
        """


def get_tab_style():
    """Returns the CSS style for the tabs."""
    return """
        :host(.bk-Tabs) {
            background-color: #15191c;
        }

        .bk-tab {
            border-right: 1px solid #303030;
        }

        :host(.bk-Tabs) .bk-header {
            color: #d0d0d0;
            border-bottom: 1px solid #d0d0d0;
        }

        .bk-tab.bk-active {
            background-color: #d0d0d0;
            color: black;
        }

        .bk-tab:not(.bk-active):hover {
            background-color: #303030;
        }

        .bk-tab:focus {
            outline: none;
        }
        """


def get_global_style():
    """Returns the global CSS style."""
    return """
        body { background: #15191c; }
        div { max-height: 96vh; max-width: 100vw; }
        """

# def get_hover_code():
#     return """
#         console.log('Current state before toggle:', sessionStorage.getItem('both_active'));
#         var current_state = sessionStorage.getItem('both_active') || 'true';
#         var new_state = current_state === 'true' ? 'false' : 'true';
#         sessionStorage.setItem('both_active', new_state);
#         console.log('New state after toggle:', sessionStorage.getItem('both_active'));

#         var toggle_to_advanced = current_state === 'true';
#         console.log('Toggle to advanced:', toggle_to_advanced);

#         function updateHoverDescriptions(hovers, normalDescs, advancedDescs) {
#             hovers.forEach(function(hover, index) {
#                 if (index < normalDescs.length && index < advancedDescs.length) {
#                     console.log('Updating hover', index);
#                     hover.tooltips = toggle_to_advanced ? advancedDescs[index] : normalDescs[index];
#                 }
#             });
#         }

#         updateHoverDescriptions(plot1_hovers, hover1_desc1, hover1_desc2);
#         updateHoverDescriptions(plot2_hovers, hover2_desc1, hover2_desc2);
#         console.log(plot1_hovers);
#         console.log(plot2_hovers);
#     """

def get_hover_code():
    return """
        var current_state = sessionStorage.getItem('advanced_mode_active') || 'false';
        var desc1, desc2;
        if (current_state === 'false') {
            desc1 = hover1_desc2;
            desc2 = hover2_desc2;
            sessionStorage.setItem('advanced_mode_active', 'true');
        } else {
            desc1 = hover1_desc1;
            desc2 = hover2_desc1;
            sessionStorage.setItem('advanced_mode_active', 'false');
        }

        for (var i = 0; i < plot1_hovers.length; i++) {
            plot1_hovers[i].tooltips = desc1[i];
        }
        for (var i = 0; i < plot2_hovers.length; i++) {
            plot2_hovers[i].tooltips = desc2[i];
        }
        """
