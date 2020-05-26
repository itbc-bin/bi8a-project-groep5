$(document).ready(function () {
    // dropdown functionality
    $('.label.ui.dropdown').dropdown();

    // table sortable
    $('.table.ui.sortable').tablesort();

    // request for co occurrence algorithm
    $(".test-button").on("click", () => {
        const id = '_' + Math.random().toString(36).substr(2, 9);
        const options = $('#multiple_select option:selected').toArray().map(item => item.text);
        const jobTitle = $('#job-title').val();
        const sendMail = $('#notify').is(":checked");
        const data = {options: options, jobTitle: jobTitle, sendMail: sendMail};
        console.log(data)
        $.ajax({
            type: 'POST',
            dataType: 'json',
            data: {data: JSON.stringify(data)},
            url: `/results/${id}`,
            async: true,
            success: [function (response) {
                $(".extra-link").remove();
                $(".form").append(`<p class="extra-link">URL: <a href="${response.url}" target="_blank">Results</a></p>`);
            }],
            error: function (response) {
                $(".form").append(`<p>something went wrong 😞</p>`)
                console.error(response)
            }
        })
    });
});


$('.trigger.example .accordion')
    .accordion({
        selector: {
            trigger: '.title .icon'
        }
    });

$(function () {
    $('#rangestart').calendar({type: 'date'});
});


// Sortable table (Van 28 tot 161)

/*
	A simple, lightweight jQuery plugin for creating sortable tables.
	https://github.com/kylefox/jquery-tablesort
	Version 0.0.11
*/

(function ($) {
    $.tablesort = function ($table, settings) {
        const self = this;
        this.$table = $table;
        this.$thead = this.$table.find('thead');
        this.settings = $.extend({}, $.tablesort.defaults, settings);
        this.$sortCells = this.$thead.length > 0 ? this.$thead.find('th:not(.no-sort)') : this.$table.find('th:not(.no-sort)');
        this.$sortCells.on('click.tablesort', function () {
            self.sort($(this));
        });
        this.index = null;
        this.$th = null;
        this.direction = null;
    };

    $.tablesort.prototype = {

        sort: function (th, direction) {
            const start = new Date(),
                self = this,
                table = this.$table,
                rowsContainer = table.find('tbody').length > 0 ? table.find('tbody') : table,
                rows = rowsContainer.find('tr').has('td, th'),
                cells = rows.find(':nth-child(' + (th.index() + 1) + ')').filter('td, th'),
                sortBy = th.data().sortBy,
                sortedMap = [];

            const unsortedValues = cells.map(function (idx, cell) {
                if (sortBy)
                    return (typeof sortBy === 'function') ? sortBy($(th), $(cell), self) : sortBy;
                return ($(this).data().sortValue != null ? $(this).data().sortValue : $(this).text());
            });
            if (unsortedValues.length === 0) return;

            //click on a different column
            if (this.index !== th.index()) {
                this.direction = 'asc';
                this.index = th.index();
            } else if (direction !== 'asc' && direction !== 'desc')
                this.direction = this.direction === 'asc' ? 'desc' : 'asc';
            else
                this.direction = direction;

            direction = this.direction === 'asc' ? 1 : -1;

            self.$table.trigger('tablesort:start', [self]);
            self.log("Sorting by " + this.index + ' ' + this.direction);

            // Try to force a browser redraw
            self.$table.css("display");
            // Run sorting asynchronously on a timeout to force browser redraw after
            // `tablesort:start` callback. Also avoids locking up the browser too much.
            setTimeout(function () {
                self.$sortCells.removeClass(self.settings.asc + ' ' + self.settings.desc);
                for (let i = 0, length = unsortedValues.length; i < length; i++) {
                    sortedMap.push({
                        index: i,
                        cell: cells[i],
                        row: rows[i],
                        value: unsortedValues[i]
                    });
                }

                sortedMap.sort(function (a, b) {
                    return self.settings.compare(a.value, b.value) * direction;
                });

                $.each(sortedMap, function (i, entry) {
                    rowsContainer.append(entry.row);
                });

                th.addClass(self.settings[self.direction]);

                self.log('Sort finished in ' + ((new Date()).getTime() - start.getTime()) + 'ms');
                self.$table.trigger('tablesort:complete', [self]);
                //Try to force a browser redraw
                self.$table.css("display");
            }, unsortedValues.length > 2000 ? 200 : 10);
        },

        log: function (msg) {
            if (($.tablesort.DEBUG || this.settings.debug) && console && console.log) {
                console.log('[tablesort] ' + msg);
            }
        },

        destroy: function () {
            this.$sortCells.off('click.tablesort');
            this.$table.data('tablesort', null);
            return null;
        }

    };

    $.tablesort.DEBUG = false;

    $.tablesort.defaults = {
        debug: $.tablesort.DEBUG,
        asc: 'sorted ascending',
        desc: 'sorted descending',
        compare: function (a, b) {
            if (a > b) {
                return 1;
            } else if (a < b) {
                return -1;
            } else {
                return 0;
            }
        }
    };

    $.fn.tablesort = function (settings) {
        let table, sortable, previous;
        return this.each(function () {
            table = $(this);
            previous = table.data('tablesort');
            if (previous) {
                previous.destroy();
            }
            table.data('tablesort', new $.tablesort(table, settings));
        });
    };

})(window.Zepto || window.jQuery);


$(function () {
    $('#submit').click(() => {
        const form_data = new FormData($('#upload-form')[0]);
        $.ajax({
            type: 'POST',
            url: '/upload_file',
            data: form_data,
            contentType: false,
            processData: false,
            dataType: 'json'
        }).done((data) => {
            alert("file succesfully uploaded!")
        })
    });
});


// sort function written by Yaris

// const tableHeaders = document.querySelectorAll(".content-table th")
// tableHeaders.forEach(headerCell => {
//     headerCell.addEventListener("click", () => {
//         const tableElement = headerCell.parentElement.parentElement.parentElement;
//         const headerIndex = Array.prototype.indexOf.call(headerCell.parentElement.children, headerCell);
//         const currentIsAscending = headerCell.classList.contains("th-sort-asc");
//
//         sortTableByColumn(tableElement, headerIndex, !currentIsAscending);
//     });
// });
//
//
// function sortTableByColumn(table, column, asc = true) {
//     const dirModifier = asc ? 1 : -1;
//     const tBody = table.tBodies[0];
//     const rows = Array.from(tBody.querySelectorAll("tr"));
//
//     // Sort each row
//     const sortedRows = rows.sort((a, b) => {
//         const aColText = a.querySelector(`td:nth-child(${column + 1})`).textContent.trim();
//         const bColText = b.querySelector(`td:nth-child(${column + 1})`).textContent.trim();
//
//         return aColText > bColText ? (dirModifier) : (-1 * dirModifier);
//     });
//
//     // Remove all existing TRs from the table
//     while (tBody.firstChild) {
//         tBody.removeChild(tBody.firstChild);
//     }
//
//     // Re-add the newly sorted rows
//     tBody.append(...sortedRows);
//
//     // Remember how the column is currently sorted
//     table.querySelectorAll("th").forEach(th => th.classList.remove("th-sort-asc", "th-sort-desc"));
//     table.querySelector(`th:nth-child(${column + 1})`).classList.toggle("th-sort-asc", asc);
//     table.querySelector(`th:nth-child(${column + 1})`).classList.toggle("th-sort-desc", !asc);
// }