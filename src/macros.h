//
// Created by Tobias Winkler on 26.10.21.
//

#ifndef PRODIGY_MACROS_H
#define PRODIGY_MACROS_H

/*
 * macro used exclusively for validating inputs to the public interface
 */
#define EXPECTS(pre, msg) if(!(pre)) throw std::invalid_argument(msg)

#endif //PRODIGY_MACROS_H
